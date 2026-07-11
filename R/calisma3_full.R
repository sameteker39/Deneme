gerekli <- c("dplyr", "tidyr", "ggplot2", "readr", "mirt", "mixtools", "patchwork")
eksik <- setdiff(gerekli, rownames(installed.packages()))
if (length(eksik)) install.packages(eksik, repos = "https://cloud.r-project.org")
invisible(lapply(gerekli, function(p) suppressMessages(library(p, character.only = TRUE))))

set.seed(20240601)
options(warn = 1)

IN_DIR  <- getwd()
OUT_DIR <- file.path(getwd(), "outputs")
FIG_DIR <- file.path(getwd(), "figs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

MAX_TIME  <- 1e5
GECERLI   <- 1:4
NT_MINN   <- 30
WINSOR_SD <- 2
NCHANCE   <- 4

bul <- function(ad) {
  aday <- c(file.path(IN_DIR, "data", ad), file.path(IN_DIR, ad))
  var  <- aday[file.exists(aday)]
  if (!length(var)) stop(sprintf("'%s' bulunamadi: %s", ad, paste(aday, collapse = " | ")))
  var[1]
}

e <- new.env(); load(bul("CALISMA2.Rdata"), envir = e)
d <- get(ls(e)[1], envir = e)
keys  <- readr::read_csv(bul("recovered_keys.csv"), show_col_types = FALSE)
items <- keys$madde_id
corr  <- setNames(keys$dogru_kod, keys$madde_id)
stopifnot(length(items) == length(unique(items)))
stopifnot(all(items %in% names(d)), all(paste0(items, "_R") %in% names(d)))

zl <- function(x) {
  if (inherits(x, "haven_labelled") || inherits(x, "labelled")) x <- unclass(x)
  suppressWarnings(as.numeric(x))
}
id_all <- as.character(zl(d$IDSTUD))
acc_df <- data.frame(IDSTUD = id_all, stringsAsFactors = FALSE)
for (it in items) { r <- zl(d[[it]])
  acc_df[[it]] <- ifelse(r %in% GECERLI, as.integer(r == corr[[it]]), NA_integer_) }
fat_df <- data.frame(IDSTUD = id_all, stringsAsFactors = FALSE)
for (it in items) { t <- zl(d[[paste0(it, "_R")]])
  fat_df[[it]] <- ifelse(is.finite(t) & t > 0 & t < MAX_TIME, t, NA_real_) }

long <- acc_df |> tidyr::pivot_longer(-IDSTUD, names_to = "items", values_to = "acc") |>
  dplyr::left_join(
    tidyr::pivot_longer(fat_df, -IDSTUD, names_to = "items", values_to = "fat"),
    by = c("IDSTUD", "items"))

winsor_top <- function(x, k = WINSOR_SD) { ok <- is.finite(x) & x > 0
  m <- mean(x[ok]); s <- sd(x[ok]); cap <- m + k * s; ifelse(ok & x > cap, cap, x) }
nt_thr <- function(x, pct = .10, maxsec = 10) { x <- x[is.finite(x) & x > 0]
  if (length(x) < NT_MINN) return(NA_real_); min(pct * mean(x), maxsec) }

long$fatx <- ave(long$fat, long$items, FUN = function(x) winsor_top(x))
thr_tab <- long |> dplyr::filter(is.finite(fatx), fatx > 0) |>
  dplyr::group_by(items) |> dplyr::summarise(thr = nt_thr(fatx), .groups = "drop")
long <- long |> dplyr::left_join(thr_tab, by = "items") |>
  dplyr::mutate(rg = as.integer(!is.na(thr) & is.finite(fat) & fat > 0 & fat <= thr))
eng <- long |> dplyr::filter(is.finite(fat), fat > 0)

val <- eng |> dplyr::filter(!is.na(acc)) |> dplyr::group_by(rg) |>
  dplyr::summarise(n = dplyr::n(), dogruluk = round(mean(acc), 3), .groups = "drop") |>
  dplyr::mutate(durum = ifelse(rg == 1, "Rapid guess", "Solution"))
readr::write_csv(val, file.path(OUT_DIR, "T0_rg_gecerlik.csv"))
cat(sprintf("RG orani: %.4f | RG dogruluk: %.3f vs caba: %.3f\n",
            mean(eng$rg), val$dogruluk[val$rg == 1], val$dogruluk[val$rg == 0]))

acc_wide <- acc_df; rownames(acc_wide) <- acc_wide$IDSTUD; acc_wide$IDSTUD <- NULL
acc_wide <- acc_wide[, items, drop = FALSE]
rg_wide <- long |> dplyr::select(IDSTUD, items, rg) |>
  tidyr::pivot_wider(names_from = items, values_from = rg)
rg_wide <- as.data.frame(rg_wide); rownames(rg_wide) <- rg_wide$IDSTUD
rg_wide <- rg_wide[rownames(acc_wide), items, drop = FALSE]
emirt_wide <- acc_wide
emirt_wide[which(rg_wide == 1, arr.ind = TRUE)] <- NA_integer_

node <- eng |> dplyr::mutate(n1 = rg, n2 = ifelse(rg == 0, acc, NA_integer_))
n1w <- node |> dplyr::select(IDSTUD, items, n1) |>
  tidyr::pivot_wider(names_from = items, values_from = n1, names_glue = "{items}_n1")
n2w <- node |> dplyr::select(IDSTUD, items, n2) |>
  tidyr::pivot_wider(names_from = items, values_from = n2, names_glue = "{items}_n2")
it_in <- items[items %in% unique(node$items)]
irt_resp <- as.data.frame(dplyr::left_join(n1w, n2w, by = "IDSTUD"))
rownames(irt_resp) <- irt_resp$IDSTUD; irt_resp$IDSTUD <- NULL
irt_resp <- irt_resp[, c(paste0(it_in, "_n1"), paste0(it_in, "_n2")), drop = FALSE]

fit_std <- mirt::mirt(acc_wide, 1, itemtype = "2PL", verbose = FALSE)
fit_em  <- mirt::mirt(emirt_wide, 1, itemtype = "2PL", verbose = FALSE)
J <- length(it_in)
spec <- mirt::mirt.model(sprintf("RG = 1-%d\n ACC = %d-%d\n COV = RG*ACC", J, J + 1, 2 * J))
fit_tree <- mirt::mirt(irt_resp, spec, itemtype = "2PL", method = "EM",
                       technical = list(NCYCLES = 2000), verbose = FALSE)

par_std <- mirt::coef(fit_std, simplify = TRUE, IRTpars = TRUE)$items
par_em  <- mirt::coef(fit_em,  simplify = TRUE, IRTpars = TRUE)$items
std_ab  <- data.frame(item = rownames(par_std), a_std = par_std[, "a"], b_std = par_std[, "b"])
em_ab   <- data.frame(item = rownames(par_em),  a_em  = par_em[, "a"],  b_em  = par_em[, "b"])
co_tree <- mirt::coef(fit_tree, simplify = TRUE)$items
acc_rows <- grepl("_n2$", rownames(co_tree))
tree_ab <- data.frame(item = sub("_n2$", "", rownames(co_tree)[acc_rows]),
                      a_tree = co_tree[acc_rows, "a2"],
                      b_tree = -co_tree[acc_rows, "d"] / co_tree[acc_rows, "a2"])
item_par <- std_ab |> dplyr::full_join(em_ab, by = "item") |>
  dplyr::full_join(tree_ab, by = "item")
readr::write_csv(item_par, file.path(OUT_DIR, "T1_madde_parametreleri.csv"))

th_std <- mirt::fscores(fit_std, method = "EAP")[, 1]
th_em  <- mirt::fscores(fit_em,  method = "EAP")[, 1]
th_tree <- mirt::fscores(fit_tree, method = "EAP")[, "ACC"]
theta <- data.frame(IDSTUD = rownames(acc_wide), theta_std = th_std, theta_em = th_em) |>
  dplyr::left_join(data.frame(IDSTUD = rownames(irt_resp), theta_tree = th_tree), by = "IDSTUD")
readr::write_csv(theta, file.path(OUT_DIR, "T2_yetenek_theta.csv"))

karsilastir <- function(x, y, etiket) {
  ok <- is.finite(x) & is.finite(y)
  data.frame(karsilastirma = etiket, n = sum(ok),
    r_pearson = round(cor(x[ok], y[ok]), 3),
    rho_spearman = round(cor(x[ok], y[ok], method = "spearman"), 3),
    ort_fark = round(mean(y[ok] - x[ok]), 3),
    RMSD = round(sqrt(mean((y[ok] - x[ok])^2)), 3))
}
cmp_all <- dplyr::bind_rows(
  karsilastir(item_par$b_std, item_par$b_em,   "b: Standart -> EM-IRT"),
  karsilastir(item_par$b_std, item_par$b_tree, "b: Standart -> IRTree(ACC)"),
  karsilastir(item_par$b_em,  item_par$b_tree, "b: EM-IRT -> IRTree(ACC)"),
  karsilastir(item_par$a_std, item_par$a_em,   "a: Standart -> EM-IRT"),
  karsilastir(item_par$a_std, item_par$a_tree, "a: Standart -> IRTree(ACC)"),
  karsilastir(item_par$a_em,  item_par$a_tree, "a: EM-IRT -> IRTree(ACC)"),
  karsilastir(theta$theta_std, theta$theta_em,   "theta: Standart -> EM-IRT"),
  karsilastir(theta$theta_std, theta$theta_tree, "theta: Standart -> IRTree(ACC)"),
  karsilastir(theta$theta_em,  theta$theta_tree, "theta: EM-IRT -> IRTree(ACC)"))
readr::write_csv(cmp_all, file.path(OUT_DIR, "T3_karsilastirma_metrikleri.csv"))
print(cmp_all)

mrtq_thr <- function(x, min_n = NT_MINN, maxrestarts = 80) {
  x <- x[is.finite(x) & x > 0]; if (length(x) < min_n) return(NA_real_)
  lx <- log(x); if (length(unique(lx)) < 2) return(NA_real_)
  out <- tryCatch(mixtools::normalmixEM(lx, k = 2, maxrestarts = maxrestarts, verb = FALSE),
                  error = function(e) NULL)
  if (is.null(out)) return(NA_real_)
  w <- out$lambda[which.min(out$mu)]
  if (w < 0.03 || w > 0.97) return(NA_real_)
  exp(stats::quantile(lx, probs = w, names = FALSE))
}
apply_rule <- function(L, method = "NT", pct = .10, maxsec = 10, fixed = NULL, winsor = TRUE) {
  L$fx <- if (winsor) ave(L$fat, L$items, FUN = function(x) winsor_top(x)) else L$fat
  th <- L |> dplyr::filter(is.finite(fx), fx > 0) |> dplyr::group_by(items) |>
    dplyr::summarise(thr = if (method == "MRTQ") mrtq_thr(fx)
                           else if (!is.null(fixed)) fixed
                           else nt_thr(fx, pct, maxsec), .groups = "drop")
  L |> dplyr::select(-dplyr::any_of("thr")) |> dplyr::left_join(th, by = "items") |>
    dplyr::mutate(rg = as.integer(!is.na(thr) & is.finite(fat) & fat > 0 & fat <= thr))
}
fit_rule <- function(Lr) {
  ew <- acc_wide
  rgw <- Lr |> dplyr::select(IDSTUD, items, rg) |>
    tidyr::pivot_wider(names_from = items, values_from = rg) |> as.data.frame()
  rownames(rgw) <- rgw$IDSTUD; rgw <- rgw[rownames(acc_wide), items, drop = FALSE]
  ew[which(rgw == 1, arr.ind = TRUE)] <- NA_integer_
  nd <- Lr |> dplyr::filter(is.finite(fat), fat > 0) |>
    dplyr::mutate(n1 = rg, n2 = ifelse(rg == 0, acc, NA_integer_))
  a1 <- nd |> dplyr::select(IDSTUD, items, n1) |>
    tidyr::pivot_wider(names_from = items, values_from = n1, names_glue = "{items}_n1")
  a2 <- nd |> dplyr::select(IDSTUD, items, n2) |>
    tidyr::pivot_wider(names_from = items, values_from = n2, names_glue = "{items}_n2")
  iin <- items[items %in% unique(nd$items)]
  ir <- as.data.frame(dplyr::left_join(a1, a2, by = "IDSTUD"))
  rownames(ir) <- ir$IDSTUD; ir$IDSTUD <- NULL
  ir <- ir[, c(paste0(iin, "_n1"), paste0(iin, "_n2")), drop = FALSE]
  Jr <- length(iin)
  sp <- mirt::mirt.model(sprintf("RG = 1-%d\n ACC = %d-%d\n COV = RG*ACC", Jr, Jr + 1, 2 * Jr))
  fe <- mirt::mirt(ew, 1, itemtype = "2PL", verbose = FALSE)
  ft <- mirt::mirt(ir, sp, itemtype = "2PL", method = "EM",
                   technical = list(NCYCLES = 2000), verbose = FALSE)
  pe <- mirt::coef(fe, simplify = TRUE, IRTpars = TRUE)$items
  co <- mirt::coef(ft, simplify = TRUE)$items; ar <- grepl("_n2$", rownames(co))
  cv <- mirt::coef(ft, simplify = TRUE)$cov
  eng2 <- Lr |> dplyr::filter(is.finite(fat), fat > 0, !is.na(acc))
  list(a_em = pe[, "a"], b_em = pe[, "b"], theta_em = mirt::fscores(fe, method = "EAP")[, 1],
       rg_rate = mean(Lr$rg[is.finite(Lr$fat) & Lr$fat > 0]),
       acc_rg = mean(eng2$acc[eng2$rg == 1]), acc_sol = mean(eng2$acc[eng2$rg == 0]),
       a_acc_tree = mean(co[ar, "a2"]),
       rho = cv["ACC", "RG"] / sqrt(cv["ACC", "ACC"] * cv["RG", "RG"]))
}
kurallar <- list(
  "NT10 (budamali, ANA)" = list(method = "NT", pct = .10, maxsec = 10, winsor = TRUE),
  "NT10 (budasiz)"       = list(method = "NT", pct = .10, maxsec = 10, winsor = FALSE),
  "NT5 (budamali)"       = list(method = "NT", pct = .05, maxsec = 10, winsor = TRUE),
  "Sabit 5 sn"           = list(method = "NT", fixed = 5, winsor = FALSE),
  "MRTQ"                 = list(method = "MRTQ", winsor = TRUE))
Lbase <- long |> dplyr::select(IDSTUD, items, acc, fat)
main_fit <- NULL; rob <- list()
for (nm in names(kurallar)) {
  fr <- do.call(apply_rule, c(list(L = Lbase), kurallar[[nm]])) |> fit_rule()
  if (is.null(main_fit)) main_fit <- fr
  rob[[nm]] <- data.frame(
    kural = nm, RG_orani = round(fr$rg_rate, 3),
    RG_dogruluk = round(fr$acc_rg, 3), caba_dogruluk = round(fr$acc_sol, 3),
    a_ort_EMIRT = round(mean(fr$a_em), 3), a_ort_IRTreeACC = round(fr$a_acc_tree, 3),
    rho_RG_ACC = round(fr$rho, 3),
    a_kor_ANA = round(cor(main_fit$a_em, fr$a_em, use = "complete.obs"), 3),
    b_kor_ANA = round(cor(main_fit$b_em, fr$b_em, use = "complete.obs"), 3),
    theta_kor_ANA = round(cor(main_fit$theta_em, fr$theta_em, use = "complete.obs"), 3))
}
rob_tab <- dplyr::bind_rows(rob)
readr::write_csv(rob_tab, file.path(OUT_DIR, "T5_duyarlilik.csv"))
print(rob_tab)

lab_map <- c("NT10 (budamali, ANA)" = "NT10\n(wins.)", "NT10 (budasiz)" = "NT10\n(unwins.)",
             "NT5 (budamali)" = "NT5", "Sabit 5 sn" = "Fixed 5 s", "MRTQ" = "MRTQ")
rob_tab$kural2 <- factor(lab_map[rob_tab$kural], levels = unname(lab_map))
POINT <- "#2a78d6"; IDENT <- "#9a9a9a"
C_RG <- "#eda100"; C_SOL <- "#1baf7a"; C_RHO <- "#4a3aa7"; C_CHANCE <- "#d1495b"
theme_apa <- function(base = 11) ggplot2::theme_classic(base_size = base) +
  ggplot2::theme(axis.title = ggplot2::element_text(face = "bold"),
    axis.line = ggplot2::element_line(color = "#444444", linewidth = 0.4),
    axis.ticks = ggplot2::element_line(color = "#444444", linewidth = 0.4),
    axis.text = ggplot2::element_text(color = "#333333"),
    legend.position = "top", legend.title = ggplot2::element_blank(),
    legend.key.size = ggplot2::unit(10, "pt"),
    plot.tag = ggplot2::element_text(face = "bold", size = 13),
    plot.margin = ggplot2::margin(6, 8, 4, 4))
sc <- function(df, xx, yy, xl, yl) {
  r <- stats::cor(df[[xx]], df[[yy]], use = "complete.obs")
  rng <- range(c(df[[xx]], df[[yy]]), na.rm = TRUE); pad <- diff(rng) * 0.06
  lims <- rng + c(-pad, pad)
  ggplot2::ggplot(df, ggplot2::aes(.data[[xx]], .data[[yy]])) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = IDENT, linewidth = 0.4) +
    ggplot2::geom_point(shape = 21, fill = POINT, color = "white", size = 2.4, stroke = 0.3, alpha = 0.9) +
    ggplot2::annotate("text", x = lims[1] + diff(lims) * 0.05, y = lims[2] - diff(lims) * 0.04,
                      label = sprintf("italic(r) == '%.3f'", r), parse = TRUE, hjust = 0, size = 3.4) +
    ggplot2::coord_equal(xlim = lims, ylim = lims) + ggplot2::labs(x = xl, y = yl) + theme_apa()
}
f1 <- (sc(item_par, "b_std", "b_em", "Difficulty (Standard 2PL)", "Difficulty (EM-IRT)") |
       sc(item_par, "b_em", "b_tree", "Difficulty (EM-IRT)", "Difficulty (IRTree)")) /
      (sc(item_par, "a_std", "a_em", "Discrimination (Standard 2PL)", "Discrimination (EM-IRT)") |
       sc(item_par, "a_em", "a_tree", "Discrimination (EM-IRT)", "Discrimination (IRTree)")) +
      patchwork::plot_annotation(tag_levels = "a")
ggplot2::ggsave(file.path(FIG_DIR, "Figure1_item_agreement.png"), f1, width = 6.8, height = 6.8, dpi = 600, bg = "white")
ggplot2::ggsave(file.path(FIG_DIR, "Figure1_item_agreement.tiff"), f1, width = 6.8, height = 6.8, dpi = 600, bg = "white", compression = "lzw")

p2a <- ggplot2::ggplot(rob_tab, ggplot2::aes(kural2, RG_orani)) +
  ggplot2::geom_col(fill = POINT, width = 0.62) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", RG_orani)), vjust = -0.5, size = 2.9, color = "#222222") +
  ggplot2::scale_y_continuous(limits = c(0, 0.17), expand = c(0, 0)) +
  ggplot2::labs(x = NULL, y = "Rapid-guessing rate") + theme_apa()
acc_long <- rob_tab |> dplyr::select(kural2, `Rapid guess` = RG_dogruluk, Solution = caba_dogruluk) |>
  tidyr::pivot_longer(-kural2, names_to = "type", values_to = "p") |>
  dplyr::mutate(type = factor(type, levels = c("Rapid guess", "Solution")))
p2b <- ggplot2::ggplot(acc_long, ggplot2::aes(kural2, p, fill = type)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(0.68), width = 0.62) +
  ggplot2::geom_hline(yintercept = 0.25, linetype = "dashed", color = C_CHANCE, linewidth = 0.5) +
  ggplot2::annotate("text", x = 4.3, y = 0.275, label = "chance = .25", color = C_CHANCE, size = 2.9, hjust = 0) +
  ggplot2::scale_fill_manual(values = c("Rapid guess" = C_RG, "Solution" = C_SOL)) +
  ggplot2::scale_y_continuous(limits = c(0, 0.65), expand = c(0, 0)) +
  ggplot2::labs(x = NULL, y = "Proportion correct") + theme_apa()
f2 <- (p2a | p2b) + patchwork::plot_annotation(tag_levels = "a")
ggplot2::ggsave(file.path(FIG_DIR, "Figure2_prevalence_validity.png"), f2, width = 7.6, height = 3.6, dpi = 600, bg = "white")
ggplot2::ggsave(file.path(FIG_DIR, "Figure2_prevalence_validity.tiff"), f2, width = 7.6, height = 3.6, dpi = 600, bg = "white", compression = "lzw")

stab <- rob_tab |> dplyr::select(kural2, Discrimination = a_kor_ANA, Difficulty = b_kor_ANA, Ability = theta_kor_ANA) |>
  tidyr::pivot_longer(-kural2, names_to = "param", values_to = "r") |>
  dplyr::mutate(param = factor(param, levels = c("Discrimination", "Difficulty", "Ability")))
p3a <- ggplot2::ggplot(stab, ggplot2::aes(kural2, r, color = param, shape = param, group = param)) +
  ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 2.6, fill = "white", stroke = 0.5) +
  ggplot2::scale_color_manual(values = c(Discrimination = "#2a78d6", Difficulty = "#1baf7a", Ability = "#eda100")) +
  ggplot2::scale_shape_manual(values = c(21, 22, 24)) +
  ggplot2::scale_y_continuous(limits = c(0.90, 1.005)) +
  ggplot2::labs(x = NULL, y = "Correlation with primary NT10") + theme_apa()
p3b <- ggplot2::ggplot(rob_tab, ggplot2::aes(kural2, rho_RG_ACC)) +
  ggplot2::geom_col(fill = C_RHO, width = 0.62) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho_RG_ACC)), vjust = 1.25, size = 2.9, color = "#222222") +
  ggplot2::geom_hline(yintercept = 0, color = "#444444", linewidth = 0.4) +
  ggplot2::scale_y_continuous(limits = c(-0.52, 0.03)) +
  ggplot2::labs(x = NULL, y = "Disengagement–ability correlation") + theme_apa()
f3 <- (p3a | p3b) + patchwork::plot_annotation(tag_levels = "a")
ggplot2::ggsave(file.path(FIG_DIR, "Figure3_sensitivity.png"), f3, width = 7.6, height = 3.6, dpi = 600, bg = "white")
ggplot2::ggsave(file.path(FIG_DIR, "Figure3_sensitivity.tiff"), f3, width = 7.6, height = 3.6, dpi = 600, bg = "white", compression = "lzw")

saveRDS(list(item_par = item_par, theta = theta, cmp = cmp_all, val = val,
             thr = thr_tab, robustness = rob_tab,
             fits = list(std = fit_std, em = fit_em, tree = fit_tree)),
        file.path(OUT_DIR, "calisma3_TUM_SONUCLAR.rds"))
cat("Bitti. Tablolar:", OUT_DIR, "| Figurler:", FIG_DIR, "\n")
