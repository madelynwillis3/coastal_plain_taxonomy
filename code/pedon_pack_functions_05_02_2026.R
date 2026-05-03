# =============================================================================
# pedon_pack_functions.R
#
# Shared, dataset-agnostic functions for building pedon_pack_150.
# Both NCSS and Perry sources call these functions after producing a
# horizon-level table that conforms to PEDON_PACK_HZ_SCHEMA below.
#
# Source this file from both the NCSS and Perry cleaning scripts.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

# -----------------------------------------------------------------------------
# REQUIRED HORIZON-TABLE SCHEMA (the "contract")
# -----------------------------------------------------------------------------
# Both NCSS and Perry must produce a horizon-level data frame containing
# AT LEAST these columns. Extra columns are allowed but ignored.
#
#   peiid                character   pedon id (unique per pedon)
#   hzdept               numeric     horizon top depth (cm)
#   hzdepb               numeric     horizon bottom depth (cm)
#   hzname               character   raw horizon designation (e.g. "Bt2", "Btg1")
#   hzname_clean         character   cleaned/squished hzname
#   hz_base              character   parsed master horizon code (e.g. "Bt", "E", "Btg")
#   L_master             character   single-letter master (A/E/B/C/O/R) or NA
#   L_suffix             character   subordinate suffix letters (t, g, x, h, k, v, p, ...)
#   texcl                character   texture class label (e.g. "Sand", "Loamy sand")
#   sandtotmeasured      numeric     total sand %
#   silttotmeasured      numeric     total silt %
#   claytotmeasured      numeric     total clay %
#   cec7                 numeric     CEC at pH 7 (cmol/kg) ; NA if unavailable
#   R, G, B              numeric     matrix color components in [0, 1] ; NA if unavailable
#
# AND a pedon-level table with one row per peiid containing taxonomy /
# location columns to be carried through:
#
#   peiid, taxonname_clean, taxorder, taxsuborder, taxgrtgroup,
#   taxsubgrp_mod1, taxsubgrp_mod2, taxsubgrp_mod3, taxpartsize,
#   POINT_X, POINT_Y   (POINT_X/Y optional; will be NA if absent)
# -----------------------------------------------------------------------------

PEDON_PACK_HZ_SCHEMA <- c(
  "peiid", "hzdept", "hzdepb",
  "hzname", "hzname_clean", "hz_base", "L_master", "L_suffix",
  "texcl",
  "sandtotmeasured", "silttotmeasured", "claytotmeasured", "cec7",
  "R", "G", "B"
)

PEDON_PACK_PED_SCHEMA <- c(
  "peiid",
  "taxonname_clean", "taxorder", "taxsuborder", "taxgrtgroup",
  "taxsubgrp_mod1", "taxsubgrp_mod2", "taxsubgrp_mod3", "taxpartsize"
)

# -----------------------------------------------------------------------------
# Utility helpers
# -----------------------------------------------------------------------------

`%||%` <- function(a, b) if (length(a) == 0 || is.na(a)) b else a

first_non_na <- function(x) {
  y <- x[!is.na(x) & x != ""]
  if (length(y) == 0) NA_character_ else as.character(y[1])
}

wmean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

overlap_thick <- function(top, bot, a, b) {
  pmax(0, pmin(bot, b) - pmax(top, a))
}

# Validate that a horizon table conforms to the schema.
validate_hz_schema <- function(hz_df, required = PEDON_PACK_HZ_SCHEMA,
                               source_name = "horizon table") {
  missing_cols <- setdiff(required, names(hz_df))
  if (length(missing_cols) > 0) {
    stop(
      "[", source_name, "] is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nProvided columns: ", paste(names(hz_df), collapse = ", ")
    )
  }
  invisible(TRUE)
}

validate_ped_schema <- function(ped_df, required = PEDON_PACK_PED_SCHEMA,
                                source_name = "pedon table") {
  missing_cols <- setdiff(required, names(ped_df))
  if (length(missing_cols) > 0) {
    stop(
      "[", source_name, "] is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nProvided columns: ", paste(names(ped_df), collapse = ", ")
    )
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Horizon-name parser (shared)
#
# Parses raw horizon designations like "Bt2", "Btg1", "2BCt", "A/E", "B't1"
# into:
#   hzname_clean  : whitespace-cleaned form
#   L_master      : first uppercase master letter (A/E/B/C/O/R)
#   L_suffix      : lowercase subordinate suffix letters (t, g, x, h, k, v, p, ...)
#   hz_base       : parsed master code used for prop_hzbase_* features
#                   - transitional masters (e.g. "BC", "AE") kept as 2 letters
#                   - otherwise master + suffix letters (e.g. "Bt", "Btg", "Ap")
#
# Used identically by NCSS and Perry so prop_hzbase_* live in the same space.
# -----------------------------------------------------------------------------

parse_one_side <- function(tok) {
  tok0 <- tok %>%
    str_replace_all("\\u00A0", " ") %>%
    str_squish() %>%
    str_replace_all("\\s+", "")
  
  if (is.na(tok0) || tok0 == "") {
    return(tibble(
      L_master = NA_character_,
      trans_master = NA_character_,
      L_suffix = NA_character_
    ))
  }
  
  tok_np <- str_remove_all(tok0, "'")
  rest   <- str_remove(tok_np, "^\\d+")  # strip lithologic discontinuity prefix
  
  trans_master <- if_else(str_detect(rest, "^[A-Z]{2}"),
                          str_sub(rest, 1, 2), NA_character_)
  
  if (!is.na(trans_master)) {
    L_master     <- str_sub(trans_master, 1, 1)
    rest2        <- str_sub(rest, 3)
    early_suffix <- ""
  } else {
    L_master <- if_else(str_detect(rest, "^[A-Z]"), str_sub(rest, 1, 1), NA_character_)
    rest2    <- if_else(!is.na(L_master), str_sub(rest, 2), rest)
    early_suffix <- str_extract(rest2, "^[a-z]+")
    early_suffix <- if_else(is.na(early_suffix), "", early_suffix)
    rest2    <- str_remove(rest2, "^[a-z]+")
  }
  
  rest3 <- str_remove(rest2, "^\\d+_\\d+")
  rest3 <- str_remove(rest3, "^\\d+")
  tail_suffix <- if_else(is.na(rest3) | rest3 == "", "", rest3)
  
  suffix_full <- paste0(early_suffix, tail_suffix)
  suffix_full <- if_else(suffix_full == "", NA_character_, suffix_full)
  
  tibble(
    L_master     = L_master,
    trans_master = trans_master,
    L_suffix     = suffix_full
  )
}

parse_hzname <- function(hz) {
  hz0   <- str_squish(str_replace_all(hz, "\\u00A0", " "))
  parts <- str_split(hz0, "\\s*/\\s*", simplify = FALSE)[[1]]
  left  <- parse_one_side(parts[1])
  
  tibble(
    hzname_clean = hz0,
    L_master     = left$L_master,
    L_suffix     = left$L_suffix,
    hz_base      = case_when(
      !is.na(left$trans_master) ~ left$trans_master,
      !is.na(left$L_master)     ~ paste0(left$L_master, coalesce(left$L_suffix, "")),
      TRUE ~ NA_character_
    )
  )
}

# Vectorized helper: take a character vector of raw horizon names, return a
# tibble with hzname (input), hzname_clean, L_master, L_suffix, hz_base.
parse_hznames <- function(hzname_vec) {
  unique_hz <- unique(hzname_vec)
  lookup <- tibble(hzname = unique_hz) %>%
    mutate(parsed = map(hzname, parse_hzname)) %>%
    unnest_wider(parsed)
  tibble(hzname = hzname_vec) %>%
    left_join(lookup, by = "hzname")
}

# -----------------------------------------------------------------------------
# Diagnostic horizon flags ("Bt-like", "sandy", "gleyed", etc.)
# Centralized so NCSS and Perry use identical rules.
# -----------------------------------------------------------------------------

is_bt_like_hz <- function(L_suffix, hzname, hzname_clean) {
  (!is.na(L_suffix) & str_detect(str_to_lower(L_suffix), "t")) |
    (!is.na(hzname_clean) & str_detect(hzname_clean, "^Bt")) |
    (!is.na(hzname) & str_detect(hzname, "^Bt"))
}

is_sandy_hz <- function(texcl, sand, clay,
                        sand_min = 70, clay_max = 15,
                        tex_fallback = c("Sand", "Loamy sand")) {
  primary  <- is.finite(sand) & is.finite(clay) & (sand >= sand_min) & (clay <= clay_max)
  fallback <- (!is.finite(sand) | !is.finite(clay)) & !is.na(texcl) &
    (str_to_lower(texcl) %in% str_to_lower(tex_fallback))
  primary | fallback
}

is_e_like_hz <- function(L_master, hz_base) {
  (!is.na(L_master)  & str_to_lower(L_master)  == "e") |
    (!is.na(hz_base) & str_to_lower(hz_base)   == "e")
}

is_gleyed_hz <- function(L_suffix) {
  !is.na(L_suffix) & str_detect(str_to_lower(L_suffix), "g")
}

is_plinthitic_hz <- function(L_suffix, hzname_clean) {
  has_v_suffix <- !is.na(L_suffix) & str_detect(str_to_lower(L_suffix), "v")
  # Tighter fallback: only "v" at end of cleaned hzname when no parsed suffix
  has_v_name <- is.na(L_suffix) & !is.na(hzname_clean) &
    str_detect(str_to_lower(hzname_clean), "v[0-9]*$|v[a-z]*[0-9]*$")
  has_v_suffix | has_v_name
}

# -----------------------------------------------------------------------------
# Per-pedon index calculators
# -----------------------------------------------------------------------------

# 1) Depth-to and thickness indices for Bt / argillic / E / gley
calc_hz_indices <- function(df, max_depth_cm = 150) {
  
  d <- df %>%
    mutate(
      hzdept = as.numeric(hzdept),
      hzdepb = as.numeric(hzdepb),
      top    = pmax(hzdept, 0, na.rm = TRUE),
      bot    = pmin(hzdepb, max_depth_cm, na.rm = TRUE),
      thk_window = pmax(0, bot - top),
      
      is_Bt       = is_bt_like_hz(L_suffix, hzname, hzname_clean),
      is_argillic = !is.na(L_suffix) & str_detect(str_to_lower(L_suffix), "t"),
      is_E        = is_e_like_hz(L_master, hz_base),
      is_gleyed   = is_gleyed_hz(L_suffix)
    ) %>%
    filter(is.finite(hzdept), is.finite(hzdepb), hzdepb > hzdept)
  
  if (nrow(d) == 0) {
    return(tibble(
      depth_to_Bt_cm = NA_real_, depth_to_argillic_cm = NA_real_,
      depth_to_E_cm = NA_real_, depth_to_gley_cm = NA_real_,
      bt_thk_0_150_cm = 0, argillic_thk_0_150_cm = 0,
      has_Bt_0_150 = FALSE, has_argillic_0_150 = FALSE,
      has_E_0_150 = FALSE, has_gley_0_150 = FALSE
    ))
  }
  
  tibble(
    depth_to_Bt_cm = if (any(d$is_Bt & is.finite(d$hzdept), na.rm = TRUE))
      min(d$hzdept[d$is_Bt], na.rm = TRUE) else NA_real_,
    depth_to_argillic_cm = if (any(d$is_argillic & is.finite(d$hzdept), na.rm = TRUE))
      min(d$hzdept[d$is_argillic], na.rm = TRUE) else NA_real_,
    depth_to_E_cm = if (any(d$is_E & is.finite(d$hzdept), na.rm = TRUE))
      min(d$hzdept[d$is_E], na.rm = TRUE) else NA_real_,
    depth_to_gley_cm = if (any(d$is_gleyed & is.finite(d$hzdept), na.rm = TRUE))
      min(d$hzdept[d$is_gleyed], na.rm = TRUE) else NA_real_,
    
    bt_thk_0_150_cm       = sum(d$thk_window[d$is_Bt],       na.rm = TRUE),
    argillic_thk_0_150_cm = sum(d$thk_window[d$is_argillic], na.rm = TRUE),
    
    has_Bt_0_150       = any(d$is_Bt       & d$thk_window > 0, na.rm = TRUE),
    has_argillic_0_150 = any(d$is_argillic & d$thk_window > 0, na.rm = TRUE),
    has_E_0_150        = any(d$is_E        & d$thk_window > 0, na.rm = TRUE),
    has_gley_0_150     = any(d$is_gleyed   & d$thk_window > 0, na.rm = TRUE)
  )
}

# 2) Sandy mantle / coarse-over-fine indices
calc_sandy_indices <- function(df, interval_top = 0, interval_bot = 150) {
  
  d <- df %>%
    mutate(
      hzdept   = as.numeric(hzdept),
      hzdepb   = as.numeric(hzdepb),
      hz_thick = hzdepb - hzdept,
      sand     = as.numeric(sandtotmeasured),
      clay     = as.numeric(claytotmeasured),
      sandy_hz = is_sandy_hz(texcl, sand, clay),
      bt_like  = is_bt_like_hz(L_suffix, hzname, hzname_clean)
    ) %>%
    filter(is.finite(hzdept), is.finite(hzdepb), hzdepb > hzdept) %>%
    arrange(hzdept, hzdepb)
  
  empty_out <- tibble(
    sandy_mantle_cm = NA_real_,
    depth_to_first_nonsandy_cm = NA_real_,
    sandy_layers_n_surface = NA_real_,
    pct_sandy_0_150 = NA_real_,
    mean_sand_0_150 = NA_real_,
    mean_clay_0_150 = NA_real_,
    has_bt_below_mantle = NA,
    delta_clay_at_mantle_break = NA_real_,
    delta_sand_at_mantle_break = NA_real_,
    abrupt_sand_to_bt_transition = NA
  )
  
  if (nrow(d) == 0) return(empty_out)
  
  surface <- min(d$hzdept, na.rm = TRUE)
  first_nonsandy_idx <- which(!d$sandy_hz)[1]
  
  sandy_mantle_cm <- if (is.na(first_nonsandy_idx)) {
    max(d$hzdepb, na.rm = TRUE) - surface
  } else {
    d$hzdept[first_nonsandy_idx] - surface
  }
  
  depth_to_first_nonsandy_cm <- if (is.na(first_nonsandy_idx)) NA_real_ else d$hzdept[first_nonsandy_idx]
  sandy_layers_n_surface     <- if (is.na(first_nonsandy_idx)) nrow(d) else (first_nonsandy_idx - 1)
  
  d_int <- d %>%
    mutate(
      ov       = overlap_thick(hzdept, hzdepb, interval_top, interval_bot),
      ov_sandy = ifelse(sandy_hz, ov, 0)
    ) %>%
    filter(ov > 0)
  
  int_thick <- sum(d_int$ov)
  pct_sandy <- if (int_thick > 0) 100 * sum(d_int$ov_sandy) / int_thick else NA_real_
  
  mean_sand <- wmean_safe(d_int$sand, d_int$ov)
  mean_clay <- wmean_safe(d_int$clay, d_int$ov)
  
  if (!is.na(first_nonsandy_idx) && first_nonsandy_idx > 1) {
    last_sandy <- d[first_nonsandy_idx - 1, ]
    first_non  <- d[first_nonsandy_idx, ]
    
    has_bt_below_mantle <- isTRUE(first_non$bt_like)
    delta_clay <- first_non$clay - last_sandy$clay
    delta_sand <- last_sandy$sand - first_non$sand
    
    abrupt <- isTRUE(has_bt_below_mantle) &&
      is.finite(delta_clay) && is.finite(delta_sand) &&
      delta_clay >= 15 && delta_sand >= 20
  } else {
    has_bt_below_mantle <- FALSE
    delta_clay <- NA_real_
    delta_sand <- NA_real_
    abrupt <- FALSE
  }
  
  tibble(
    sandy_mantle_cm = sandy_mantle_cm,
    depth_to_first_nonsandy_cm = depth_to_first_nonsandy_cm,
    sandy_layers_n_surface = sandy_layers_n_surface,
    pct_sandy_0_150 = pct_sandy,
    mean_sand_0_150 = mean_sand,
    mean_clay_0_150 = mean_clay,
    has_bt_below_mantle = has_bt_below_mantle,
    delta_clay_at_mantle_break = delta_clay,
    delta_sand_at_mantle_break = delta_sand,
    abrupt_sand_to_bt_transition = abrupt
  )
}

# 3) Plinthite indices
calc_plinthite_indices <- function(df, upper_limit = 150) {
  
  d <- df %>%
    mutate(
      hzdept = as.numeric(hzdept),
      hzdepb = as.numeric(hzdepb),
      hz_thick = hzdepb - hzdept,
      plinthitic_hz = is_plinthitic_hz(L_suffix, hzname_clean)
    ) %>%
    filter(is.finite(hzdept), is.finite(hzdepb), hzdepb > hzdept) %>%
    arrange(hzdept, hzdepb)
  
  empty_out <- tibble(
    plinthite_flag = NA, top_depth_plinthite = NA_real_,
    plinthite_thickness_cm = NA_real_,
    pct_plinthite_by_hz = NA_real_,
    pct_plinthite_by_thickness = NA_real_,
    plinthite_overlap_0_150_cm = NA_real_,
    plinthite_within_150cm_flag = NA
  )
  
  if (nrow(d) == 0) return(empty_out)
  
  profile_top       <- min(d$hzdept, na.rm = TRUE)
  profile_bottom    <- max(d$hzdepb, na.rm = TRUE)
  profile_thickness <- profile_bottom - profile_top
  
  d <- d %>%
    mutate(ov_0_150 = overlap_thick(hzdept, hzdepb, 0, upper_limit))
  
  plinthite_flag         <- any(d$plinthitic_hz, na.rm = TRUE)
  top_depth_plinthite    <- if (plinthite_flag) min(d$hzdept[d$plinthitic_hz], na.rm = TRUE) else NA_real_
  plinthite_thickness_cm <- sum(d$hz_thick[d$plinthitic_hz], na.rm = TRUE)
  pct_plinthite_by_hz    <- 100 * mean(d$plinthitic_hz, na.rm = TRUE)
  
  pct_plinthite_by_thickness <- if (is.finite(profile_thickness) && profile_thickness > 0) {
    100 * plinthite_thickness_cm / profile_thickness
  } else NA_real_
  
  plinthite_overlap_0_150_cm  <- sum(d$ov_0_150[d$plinthitic_hz], na.rm = TRUE)
  plinthite_within_150cm_flag <- any(d$plinthitic_hz & d$ov_0_150 > 0, na.rm = TRUE)
  
  tibble(
    plinthite_flag = plinthite_flag,
    top_depth_plinthite = top_depth_plinthite,
    plinthite_thickness_cm = plinthite_thickness_cm,
    pct_plinthite_by_hz = pct_plinthite_by_hz,
    pct_plinthite_by_thickness = pct_plinthite_by_thickness,
    plinthite_overlap_0_150_cm = plinthite_overlap_0_150_cm,
    plinthite_within_150cm_flag = plinthite_within_150cm_flag
  )
}

# 4) Argillic / Kandic / Pale-style indices
calc_argillic_kandic_pale_indices <- function(df, upper = 150,
                                              clay_jump_ratio = 1.2,
                                              clay_jump_abs = 3) {
  
  d <- df %>%
    mutate(
      hzdept   = as.numeric(hzdept),
      hzdepb   = as.numeric(hzdepb),
      hz_thick = hzdepb - hzdept,
      clay     = as.numeric(claytotmeasured),
      sand     = as.numeric(sandtotmeasured),
      cec7     = as.numeric(cec7),
      bt_like  = is_bt_like_hz(L_suffix, hzname, hzname_clean),
      ov_0_150 = overlap_thick(hzdept, hzdepb, 0, upper),
      in_0_150 = ov_0_150 > 0
    ) %>%
    filter(is.finite(hzdept), is.finite(hzdepb), hzdepb > hzdept) %>%
    arrange(hzdept, hzdepb)
  
  empty_out <- tibble(
    has_bt_like = NA, bt_top_cm = NA_real_, bt_thickness_cm = NA_real_,
    max_clay_ratio_0_150 = NA_real_, max_clay_delta_0_150 = NA_real_,
    depth_to_first_clay_jump = NA_real_,
    argillic_top_cm = NA_real_, argillic_thickness_cm = NA_real_,
    argillic_contiguous_max_cm = NA_real_,
    pct_hz_with_cec7 = NA_real_,
    min_cec_per_kg_clay_0_150 = NA_real_,
    median_cec_per_kg_clay_bt_zone = NA_real_,
    has_major_clay_drop_within_150 = NA,
    depth_to_20pct_clay_decrease = NA_real_,
    clay_persistence_index_0_150 = NA_real_
  )
  
  if (nrow(d) == 0) return(empty_out)
  
  has_bt_like     <- any(d$bt_like, na.rm = TRUE)
  bt_top_cm       <- if (has_bt_like) min(d$hzdept[d$bt_like], na.rm = TRUE) else NA_real_
  bt_thickness_cm <- sum(d$hz_thick[d$bt_like], na.rm = TRUE)
  
  d150 <- d %>% filter(in_0_150)
  
  if (!any(is.finite(d150$clay))) {
    max_clay_ratio   <- NA_real_
    max_clay_delta   <- NA_real_
    depth_first_jump <- NA_real_
    arg_top          <- NA_real_
    arg_thick        <- NA_real_
    arg_contig       <- NA_real_
  } else {
    clay_vals <- d150$clay
    tops      <- d150$hzdept
    
    run_min_above <- lag(cummin(clay_vals), default = NA_real_)
    ratio <- clay_vals / run_min_above
    delta <- clay_vals - run_min_above
    
    jump <- is.finite(ratio) & is.finite(delta) &
      (ratio >= clay_jump_ratio | delta >= clay_jump_abs)
    
    max_clay_ratio <- suppressWarnings(max(ratio[jump], na.rm = TRUE))
    if (!is.finite(max_clay_ratio)) max_clay_ratio <- NA_real_
    
    max_clay_delta <- suppressWarnings(max(delta[jump], na.rm = TRUE))
    if (!is.finite(max_clay_delta)) max_clay_delta <- NA_real_
    
    depth_first_jump <- if (any(jump, na.rm = TRUE)) tops[which(jump)[1]] else NA_real_
    
    arg_zone <- rep(FALSE, nrow(d150))
    if (is.finite(depth_first_jump)) {
      arg_zone <- d150$hzdept >= depth_first_jump & (d150$bt_like | jump)
    }
    
    arg_top   <- if (any(arg_zone)) min(d150$hzdept[arg_zone], na.rm = TRUE) else NA_real_
    arg_thick <- sum(d150$hz_thick[arg_zone], na.rm = TRUE)
    
    if (any(arg_zone)) {
      z <- d150 %>%
        mutate(is_arg = arg_zone) %>%
        mutate(run = cumsum(is_arg != lag(is_arg, default = first(is_arg))))
      arg_contig <- z %>%
        filter(is_arg) %>%
        group_by(run) %>%
        summarise(cm = sum(hz_thick, na.rm = TRUE), .groups = "drop") %>%
        summarise(mx = max(cm, na.rm = TRUE)) %>%
        pull(mx)
    } else {
      arg_contig <- NA_real_
    }
  }
  
  pct_hz_with_cec7 <- mean(is.finite(d150$cec7), na.rm = TRUE) * 100
  
  d150 <- d150 %>%
    mutate(
      cec_per_kg_clay = ifelse(
        is.finite(cec7) & is.finite(clay) & clay > 0,
        cec7 / (clay / 100),
        NA_real_
      )
    )
  
  min_cec_per_kg_clay <- suppressWarnings(min(d150$cec_per_kg_clay, na.rm = TRUE))
  if (!is.finite(min_cec_per_kg_clay)) min_cec_per_kg_clay <- NA_real_
  
  median_cec_per_kg_clay_bt <- suppressWarnings(
    median(d150$cec_per_kg_clay[d150$bt_like], na.rm = TRUE)
  )
  if (!is.finite(median_cec_per_kg_clay_bt)) median_cec_per_kg_clay_bt <- NA_real_
  
  if (any(is.finite(d150$clay))) {
    i_max        <- which.max(d150$clay)
    clay_max     <- d150$clay[i_max]
    depth_at_max <- d150$hzdept[i_max]
    
    below <- d150 %>% filter(hzdept > depth_at_max)
    
    drop20 <- if (nrow(below) > 0 && is.finite(clay_max) && clay_max > 0) {
      (clay_max - below$clay) / clay_max >= 0.20
    } else FALSE
    
    has_major_drop <- any(drop20, na.rm = TRUE)
    depth_to_drop  <- if (has_major_drop) below$hzdept[which(drop20)[1]] else NA_real_
    
    zone_for_persist <- d150 %>%
      filter(is.finite(clay)) %>%
      filter(bt_like | (!is.na(arg_top) & hzdept >= arg_top))
    
    persist <- if (nrow(zone_for_persist) >= 2 && mean(zone_for_persist$clay) > 0) {
      1 - (sd(zone_for_persist$clay) / mean(zone_for_persist$clay))
    } else NA_real_
  } else {
    has_major_drop <- NA
    depth_to_drop  <- NA_real_
    persist        <- NA_real_
  }
  
  tibble(
    has_bt_like = has_bt_like,
    bt_top_cm = bt_top_cm,
    bt_thickness_cm = bt_thickness_cm,
    max_clay_ratio_0_150 = max_clay_ratio,
    max_clay_delta_0_150 = max_clay_delta,
    depth_to_first_clay_jump = depth_first_jump,
    argillic_top_cm = arg_top,
    argillic_thickness_cm = arg_thick,
    argillic_contiguous_max_cm = arg_contig,
    pct_hz_with_cec7 = pct_hz_with_cec7,
    min_cec_per_kg_clay_0_150 = min_cec_per_kg_clay,
    median_cec_per_kg_clay_bt_zone = median_cec_per_kg_clay_bt,
    has_major_clay_drop_within_150 = has_major_drop,
    depth_to_20pct_clay_decrease = depth_to_drop,
    clay_persistence_index_0_150 = persist
  )
}

# -----------------------------------------------------------------------------
# Master assembler: hz_table + ped_table -> pedon_pack_150
# -----------------------------------------------------------------------------
#
# Inputs:
#   hz_table   : horizon-level df conforming to PEDON_PACK_HZ_SCHEMA
#   ped_table  : pedon-level df conforming to PEDON_PACK_PED_SCHEMA
#                (one row per peiid; carries taxonomy + optional location)
#   max_depth_cm : depth window (default 150)
#
# Output: pedon_pack tibble, one row per peiid
# -----------------------------------------------------------------------------
build_pedon_pack <- function(hz_table, ped_table,
                             max_depth_cm = 150,
                             source_name = "dataset") {
  
  validate_hz_schema(hz_table, source_name = paste0(source_name, " horizon table"))
  validate_ped_schema(ped_table, source_name = paste0(source_name, " pedon table"))
  
  hz <- hz_table %>%
    mutate(peiid = as.character(peiid))
  
  ped <- ped_table %>%
    mutate(peiid = as.character(peiid)) %>%
    distinct(peiid, .keep_all = TRUE)
  
  # ---- 1. Weighted means within 0–max_depth_cm window ---------------------
  pack_wmeans <- hz %>%
    mutate(
      top = pmax(as.numeric(hzdept), 0),
      bot = pmin(as.numeric(hzdepb), max_depth_cm),
      thk = pmax(0, bot - top)
    ) %>%
    filter(thk > 0) %>%
    group_by(peiid) %>%
    summarise(
      sand_w = wmean_safe(as.numeric(sandtotmeasured), thk),
      silt_w = wmean_safe(as.numeric(silttotmeasured), thk),
      clay_w = wmean_safe(as.numeric(claytotmeasured), thk),
      cec7_w = wmean_safe(as.numeric(cec7), thk),
      R_w    = wmean_safe(as.numeric(R), thk),
      G_w    = wmean_safe(as.numeric(G), thk),
      B_w    = wmean_safe(as.numeric(B), thk),
      
      any_t_suffix = any(str_detect(str_to_lower(coalesce(L_suffix, "")), "t"), na.rm = TRUE),
      any_g_suffix = any(str_detect(str_to_lower(coalesce(L_suffix, "")), "g"), na.rm = TRUE),
      any_x_suffix = any(str_detect(str_to_lower(coalesce(L_suffix, "")), "x"), na.rm = TRUE),
      any_k_suffix = any(str_detect(str_to_lower(coalesce(L_suffix, "")), "k"), na.rm = TRUE),
      
      .groups = "drop"
    )
  
  # ---- 2. Horizon-base proportions in 0-max_depth_cm window ---------------
  hz_props <- hz %>%
    mutate(
      top = pmax(as.numeric(hzdept), 0),
      bot = pmin(as.numeric(hzdepb), max_depth_cm),
      thk = pmax(0, bot - top),
      hz_base2 = if_else(is.na(hz_base) | hz_base == "", "Unknown", hz_base)
    ) %>%
    filter(thk > 0) %>%
    group_by(peiid, hz_base2) %>%
    summarise(thk = sum(thk), .groups = "drop") %>%
    group_by(peiid) %>%
    mutate(prop = thk / sum(thk)) %>%
    ungroup() %>%
    select(peiid, hz_base2, prop) %>%
    pivot_wider(
      names_from   = hz_base2,
      values_from  = prop,
      values_fill  = 0,
      names_prefix = "prop_hzbase_"
    )
  
  # ---- 3. Per-pedon index calculations ------------------------------------
  hz_indices <- hz %>%
    group_by(peiid) %>%
    group_modify(~ calc_hz_indices(.x, max_depth_cm = max_depth_cm)) %>%
    ungroup()
  
  sandy_indices <- hz %>%
    group_by(peiid) %>%
    group_modify(~ calc_sandy_indices(.x, 0, max_depth_cm)) %>%
    ungroup()
  
  plinthite_indices <- hz %>%
    group_by(peiid) %>%
    group_modify(~ calc_plinthite_indices(.x, upper_limit = max_depth_cm)) %>%
    ungroup()
  
  arg_kand_pale_indices <- hz %>%
    group_by(peiid) %>%
    group_modify(~ calc_argillic_kandic_pale_indices(.x, upper = max_depth_cm)) %>%
    ungroup()
  
  # ---- 4. Stitch everything onto the pedon-level taxonomy table -----------
  pedon_pack <- ped %>%
    left_join(pack_wmeans,           by = "peiid") %>%
    left_join(hz_props,              by = "peiid") %>%
    left_join(hz_indices,            by = "peiid") %>%
    left_join(sandy_indices,         by = "peiid") %>%
    left_join(plinthite_indices,     by = "peiid") %>%
    left_join(arg_kand_pale_indices, by = "peiid")
  
  pedon_pack
}
