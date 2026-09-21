library(ggplot2)
library(patchwork)
library(scales)

# ---- 1. Run simulations, extracting data only (no base-R plots) ----
ctrl_res <- run(phase_def = phase_def, n = 1, eta = eta, lambda = lambda,
                alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
                gamma = gamma, delta = delta, sigma0 = sigma0,
                plot = FALSE, plot_raw = FALSE)[[1]]

ptsd_res <- run(phase_def = phase_def, n = 1, eta = eta_ptsd, lambda = lambda,
                alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
                gamma = gamma, delta = delta, sigma0 = sigma0,
                plot = FALSE, plot_raw = FALSE)[[1]]

# ---- 2. Tidy-data helpers ----

# Overall (belief-weighted) V, summed across all cues, for panels A/B
make_overall_df <- function(res) {
  cue_names <- names(res$V)
  do.call(rbind, lapply(cue_names, function(cue) {
    v <- res$V[[cue]][, 1]                       # dimension D1
    data.frame(trial = seq_along(v), cue = cue, V = v)
  }))
}

# Raw V within one latent state, for panels C/D/E
# NOTE: the associative-value trace only starts on the trial AFTER the state
# is created (e.g. if the state is created on trial 11, the V line begins at
# trial 12), since V for a not-yet-created state is meaningless. The BELIEF
# trace, however, is kept for all trials, since belief in a state is
# well-defined (and typically ~0) even before the state formally exists.
make_state_df <- function(res, state) {
  cue_names <- names(res$V_raw)
  
  # Trial on which this state was first created (Lmax reaches >= state)
  creation_trial <- which(res$Lmax >= state)[1]
  start_trial <- if (is.na(creation_trial)) Inf else creation_trial + 1
  
  v_df <- do.call(rbind, lapply(cue_names, function(cue) {
    v <- res$V_raw[[cue]][, 1, state]               # dimension D1, this state
    data.frame(trial = seq_along(v), cue = cue, V = v)
  }))
  v_df <- v_df[v_df$trial >= start_trial, ]
  
  belief <- res$lsb[, state]
  belief_df <- data.frame(trial = seq_along(belief), belief = belief)
  # belief_df is NOT filtered — kept for all trials
  
  list(v_df = v_df, belief = belief_df)
}

# Phase-transition trial indices (vertical dashed reference lines)
phase_bounds <- function() {
  ends <- cumsum(sapply(phase_def, function(ph) sum(unlist(ph$trials))))
  ends[-length(ends)] + 0.5
}

# One shared color per cue so it stays consistent across all panels
cue_pal    <- c(A = "#1b9e77")   # add more entries here if you add cues
belief_col <- "grey30"

# ---- 3. Panel builders ----

plot_overall <- function(res, tag, subtitle) {
  df <- make_overall_df(res)
  ggplot(df, aes(trial, V, color = cue)) +
    geom_vline(xintercept = phase_bounds(), linetype = "dashed", color = "grey70") +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = cue_pal, guide = "none") +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    labs(subtitle = subtitle, x = "Trial", y = "Associative value (V)", tag = tag) +
    theme_classic(base_size = 11) +
    theme(plot.subtitle = element_text(hjust = 0.5, face = "italic"),
          plot.tag = element_text(face = "bold"))
}

plot_state <- function(res, state, tag, subtitle) {
  d         <- make_state_df(res, state)
  v_df      <- d$v_df
  belief_df <- d$belief
  y_range   <- c(-0.5, 0.5)
  
  # Rescale belief (0-1) onto the same axis range so it can share the panel
  belief_df$V_scaled <- y_range[1] + belief_df$belief * diff(y_range)
  
  ggplot() +
    geom_vline(xintercept = phase_bounds(), linetype = "dashed", color = "grey70") +
    geom_line(data = v_df, aes(trial, V, color = cue), linewidth = 0.9) +
    geom_line(data = belief_df, aes(trial, V_scaled),
              color = belief_col, linetype = "twodash", linewidth = 0.8) +
    scale_color_manual(values = cue_pal, guide = "none") +
    scale_x_continuous(limits = c(1, max(res$par$ntrials, 1))) +
    scale_y_continuous(
      limits = y_range,
      sec.axis = sec_axis(~ (. - y_range[1]) / diff(y_range),
                          name = "Belief in state", labels = percent)
    ) +
    labs(subtitle = subtitle, x = "Trial", y = "Associative value (V)", tag = tag) +
    theme_classic(base_size = 11) +
    theme(plot.subtitle = element_text(hjust = 0.5, face = "italic"),
          plot.tag = element_text(face = "bold"))
}

# ---- 4. Build the five panels ----
pA <- plot_overall(ctrl_res, "A", "Overall learning")
pB <- plot_overall(ptsd_res, "B", "Overall learning")
pC <- plot_state(ctrl_res, 1, "C", "Latent state 1")
pE <- plot_state(ptsd_res, 1, "E", "Latent state 1")
pD <- plot_state(ctrl_res, 2, "D", "Latent state 2")

# ---- 5. Column headers ("Control" / "PTSD") ----
header <- function(label) {
  ggplot() + theme_void() +
    annotate("text", x = 0, y = 0, label = label, fontface = "bold", size = 5.5) +
    theme(plot.margin = margin(2, 2, 2, 2))
}
hdr_ctrl <- header("Control")
hdr_ptsd <- header("PTSD")

# ---- 6. Assemble final figure ----
row_hdr <- hdr_ctrl + hdr_ptsd
row1    <- pA + pB
row2    <- pC + pE
row3    <- pD + plot_spacer()   # PTSD never forms a 2nd latent state

figure <- row_hdr / row1 / row2 / row3 +
  plot_layout(heights = c(0.1, 1, 1, 1))

# ---- 7. Save at publication resolution ----
ggsave("PTSD_vs_Control_LatentState.png", figure,
       width = 9, height = 10, dpi = 400, bg = "white")