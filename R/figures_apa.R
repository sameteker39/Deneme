# #####################################################################
# APA UYUMLU FIGURLER (Ingilizce) — ggplot2
# Girdi:  outputs/T1_madde_parametreleri.csv , outputs/T5_duyarlilik.csv
#         (calisma3_irtree_vs_emirt.R ciktilari)
# Cikti:  figs/Figure1_item_agreement.(png|tiff)
#         figs/Figure2_prevalence_validity.(png|tiff)
#         figs/Figure3_sensitivity.(png|tiff)   -- 600 dpi
#
# Tasarim: gri-tonlamaya uygun, renk-korligune guvenli palet (dogrulanmis),
# kimlik cizgisi, minimal eksen, panel harfleri (a-d). APA 7 uyumlu.
# #####################################################################

gerekli <- c("dplyr", "tidyr", "ggplot2", "readr", "patchwork")
eksik <- setdiff(gerekli, rownames(installed.packages()))
if (length(eksik)) install.packages(eksik, repos = "https://cloud.r-project.org")
invisible(lapply(gerekli, function(p) suppressMessages(library(p, character.only = TRUE))))

IN  <- file.path(getwd(), "outputs")
OUT <- file.path(getwd(), "figs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

ip <- readr::read_csv(file.path(IN, "T1_madde_parametreleri.csv"), show_col_types = FALSE)
t5 <- readr::read_csv(file.path(IN, "T5_duyarlilik.csv"), show_col_types = FALSE)

# Kural etiketlerini Ingilizce'ye cevir + sirayi sabitle
lab_map <- c("NT10 (budamali, ANA)" = "NT10\n(wins.)",
             "NT10 (budasiz)"       = "NT10\n(unwins.)",
             "NT5 (budamali)"       = "NT5",
             "Sabit 5 sn"           = "Fixed 5 s",
             "MRTQ"                 = "MRTQ")
t5$kural2 <- factor(lab_map[t5$kural], levels = unname(lab_map))

# --- Palet (dogrulanmis: CVD-guvenli) ---
POINT <- "#2a78d6"; IDENT <- "#9a9a9a"
C_RG <- "#eda100"; C_SOL <- "#1baf7a"; C_RHO <- "#4a3aa7"; C_CHANCE <- "#d1495b"

# --- APA temasi ---
theme_apa <- function(base = 11) {
  ggplot2::theme_classic(base_size = base) +
    ggplot2::theme(
      axis.title  = ggplot2::element_text(face = "bold"),
      axis.line   = ggplot2::element_line(color = "#444444", linewidth = 0.4),
      axis.ticks  = ggplot2::element_line(color = "#444444", linewidth = 0.4),
      axis.text   = ggplot2::element_text(color = "#333333"),
      legend.position = "top", legend.title = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(10, "pt"),
      plot.tag = ggplot2::element_text(face = "bold", size = 13),
      plot.margin = ggplot2::margin(6, 8, 4, 4))
}

# ===================== FIGURE 1: madde parametre ozdesligi ==========
sc <- function(df, xx, yy, xl, yl) {
  r <- stats::cor(df[[xx]], df[[yy]], use = "complete.obs")
  rng <- range(c(df[[xx]], df[[yy]]), na.rm = TRUE); pad <- diff(rng) * 0.06
  lims <- rng + c(-pad, pad)
  ggplot2::ggplot(df, ggplot2::aes(.data[[xx]], .data[[yy]])) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = IDENT, linewidth = 0.4) +
    ggplot2::geom_point(shape = 21, fill = POINT, color = "white",
                        size = 2.4, stroke = 0.3, alpha = 0.9) +
    ggplot2::annotate("text", x = lims[1] + diff(lims) * 0.05,
                      y = lims[2] - diff(lims) * 0.04,
                      label = sprintf("italic(r) == '%.3f'", r), parse = TRUE,
                      hjust = 0, size = 3.4) +
    ggplot2::coord_equal(xlim = lims, ylim = lims) +
    ggplot2::labs(x = xl, y = yl) + theme_apa()
}
f1 <- (sc(ip, "b_std", "b_em",   "Difficulty (Standard 2PL)",  "Difficulty (EM-IRT)") |
       sc(ip, "b_em",  "b_tree", "Difficulty (EM-IRT)",        "Difficulty (IRTree)")) /
      (sc(ip, "a_std", "a_em",   "Discrimination (Standard 2PL)", "Discrimination (EM-IRT)") |
       sc(ip, "a_em",  "a_tree", "Discrimination (EM-IRT)",       "Discrimination (IRTree)")) +
      patchwork::plot_annotation(tag_levels = "a")
ggplot2::ggsave(file.path(OUT, "Figure1_item_agreement.png"), f1,
                width = 6.8, height = 6.8, dpi = 600, bg = "white")
ggplot2::ggsave(file.path(OUT, "Figure1_item_agreement.tiff"), f1,
                width = 6.8, height = 6.8, dpi = 600, bg = "white", compression = "lzw")

# ===================== FIGURE 2: yayginlik + gecerlik ==============
p2a <- ggplot2::ggplot(t5, ggplot2::aes(kural2, RG_orani)) +
  ggplot2::geom_col(fill = POINT, width = 0.62) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", RG_orani)),
                     vjust = -0.5, size = 2.9, color = "#222222") +
  ggplot2::scale_y_continuous(limits = c(0, 0.17), expand = c(0, 0)) +
  ggplot2::labs(x = NULL, y = "Rapid-guessing rate") + theme_apa()

acc_long <- t5 |>
  dplyr::select(kural2, `Rapid guess` = RG_dogruluk, Solution = caba_dogruluk) |>
  tidyr::pivot_longer(-kural2, names_to = "type", values_to = "p") |>
  dplyr::mutate(type = factor(type, levels = c("Rapid guess", "Solution")))
p2b <- ggplot2::ggplot(acc_long, ggplot2::aes(kural2, p, fill = type)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(0.68), width = 0.62) +
  ggplot2::geom_hline(yintercept = 0.25, linetype = "dashed", color = C_CHANCE, linewidth = 0.5) +
  ggplot2::annotate("text", x = 4.3, y = 0.275, label = "chance = .25",
                    color = C_CHANCE, size = 2.9, hjust = 0) +
  ggplot2::scale_fill_manual(values = c("Rapid guess" = C_RG, "Solution" = C_SOL)) +
  ggplot2::scale_y_continuous(limits = c(0, 0.65), expand = c(0, 0)) +
  ggplot2::labs(x = NULL, y = "Proportion correct") + theme_apa()
f2 <- (p2a | p2b) + patchwork::plot_annotation(tag_levels = "a")
ggplot2::ggsave(file.path(OUT, "Figure2_prevalence_validity.png"), f2,
                width = 7.6, height = 3.6, dpi = 600, bg = "white")
ggplot2::ggsave(file.path(OUT, "Figure2_prevalence_validity.tiff"), f2,
                width = 7.6, height = 3.6, dpi = 600, bg = "white", compression = "lzw")

# ===================== FIGURE 3: duyarlilik ========================
stab <- t5 |>
  dplyr::select(kural2, Discrimination = a_kor_ANA, Difficulty = b_kor_ANA,
                Ability = theta_kor_ANA) |>
  tidyr::pivot_longer(-kural2, names_to = "param", values_to = "r") |>
  dplyr::mutate(param = factor(param, levels = c("Discrimination", "Difficulty", "Ability")))
p3a <- ggplot2::ggplot(stab, ggplot2::aes(kural2, r, color = param, shape = param, group = param)) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 2.6, fill = "white", stroke = 0.5) +
  ggplot2::scale_color_manual(values = c(Discrimination = "#2a78d6", Difficulty = "#1baf7a", Ability = "#eda100")) +
  ggplot2::scale_shape_manual(values = c(21, 22, 24)) +
  ggplot2::scale_y_continuous(limits = c(0.90, 1.005)) +
  ggplot2::labs(x = NULL, y = "Correlation with primary NT10") + theme_apa()

p3b <- ggplot2::ggplot(t5, ggplot2::aes(kural2, rho_RG_ACC)) +
  ggplot2::geom_col(fill = C_RHO, width = 0.62) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho_RG_ACC)),
                     vjust = 1.25, size = 2.9, color = "#222222") +
  ggplot2::geom_hline(yintercept = 0, color = "#444444", linewidth = 0.4) +
  ggplot2::scale_y_continuous(limits = c(-0.52, 0.03)) +
  ggplot2::labs(x = NULL, y = "Disengagement–ability correlation") + theme_apa()
f3 <- (p3a | p3b) + patchwork::plot_annotation(tag_levels = "a")
ggplot2::ggsave(file.path(OUT, "Figure3_sensitivity.png"), f3,
                width = 7.6, height = 3.6, dpi = 600, bg = "white")
ggplot2::ggsave(file.path(OUT, "Figure3_sensitivity.tiff"), f3,
                width = 7.6, height = 3.6, dpi = 600, bg = "white", compression = "lzw")

cat("Figurler yazildi:", OUT, "\n")
