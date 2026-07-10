# #####################################################################
# CALISMA 3 - TIMSS 2023 (Turkiye, 8. sinif matematik, coktan secmeli)
# IRTREE ile EM-IRT (Effort-Moderated IRT) MADDE ve YETENEK PARAMETRELERININ
# KARSILASTIRILMASI.  Kendi kendine yeten analiz betigi (Ctrl+A -> Run).
#
# Girdi (IN_DIR icinde):
#   - CALISMA2.Rdata        : TIMSS 2023 Turkiye genis format
#                             (ME##### ham yanit 1-4 ; ME#####_R ilk yanit suresi ;
#                              BSMMAT01-05 plausible values)
#   - recovered_keys.csv    : 43 MC madde icin ampirik cevap anahtari
#                             (sutunlar: madde_id, dogru_kod, ability_gap, n)
#
# Cikti (OUT_DIR/): parametre tablolari (.csv), karsilastirma figurleri (.png),
#   tam sonuc (.rds), konsol gunlugu (.txt).
#
# TASARIM
#   Esik      : NT10 (Wise & Ma, 2012), NT oncesi uc sureler ort+2SS'te budanir.
#   Modeller  : (1) Standart 2PL   - tum yanitlar, cabayi yok sayar (baseline)
#               (2) EM-IRT         - Wise & DeMars (2006). Hizli tahmin yanitlari
#                                    sans fonksiyonuna sabitlenir; bu, madde/yetenek
#                                    kestirimi acisindan "hizli tahmin -> NA + 2PL"
#                                    ile OZDESTIR (sabit likelihood katkisi).
#               (3) IRTree 2PL     - iki dugum: n1=RG (hizli tahmin), n2=ACC (dogruluk);
#                                    mirt cok-boyutlu, ACC dugumu effortful yanitlari
#                                    RG boyutuyla ES ZAMANLI kestirir.
#   Karsilastirma: madde gucluk (b), ayirt edicilik (a) ve yetenek (theta) uzerinden
#                  korelasyon, ortalama fark, RMSD ve sacilim grafikleri.
# #####################################################################

# ===================== 0. KURULUM ===================================
gerekli <- c("dplyr", "tidyr", "ggplot2", "mirt", "readr", "mixtools")
kurulu  <- rownames(installed.packages())
eksik   <- setdiff(gerekli, kurulu)
if (length(eksik)) install.packages(eksik, repos = "https://cloud.r-project.org")
invisible(lapply(gerekli, function(p) suppressMessages(library(p, character.only = TRUE))))

set.seed(20240601)
options(warn = 1)

IN_DIR  <- getwd()
OUT_DIR <- file.path(getwd(), "outputs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

.logcon <- file(file.path(OUT_DIR, "calisma3_konsol_log.txt"), open = "wt")
sink(.logcon, split = TRUE)
cat("##### CALISMA 3 - IRTree vs EM-IRT #####\n")
cat("Tarih:", format(Sys.time()), " | Cikti:", OUT_DIR, "\n\n")

# --- Sabitler ---
# DIKKAT: verideki eksik-sure sentinel'i 9999999999'dur; "t < 1e10" bunu DISLAMAZ
# (9999999999 < 1e10). Bu yuzden gecerli RT ust siniri 1e5 alinir (gercek RT < ~900 sn).
MAX_TIME  <- 1e5           # gecerli yanit suresi ust siniri (sentinel haric)
GECERLI   <- 1:4           # gecerli MC secenek
NT_MINN   <- 30            # esik icin madde basi en az gozlem
WINSOR_SD <- 2             # uc sure budama (ort + 2 SS)
NCHANCE   <- 4             # 4 secenekli MC -> sans = 1/4
theme_apa <- function() ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 axis.title = ggplot2::element_text(face = "bold"),
                 legend.position = "bottom", legend.title = ggplot2::element_blank())
savefig <- function(p, f, w = 6.5, h = 5)
  ggplot2::ggsave(file.path(OUT_DIR, f), p, width = w, height = h, dpi = 600, bg = "white")

# ===================== 1. VERI HAZIRLAMA ============================
cat("== 1. Veri okunuyor ve puanlaniyor ==\n")
# Dosyayi hem <klasor>/data/ altinda hem de dogrudan calisma klasorunde arar.
bul <- function(ad) {
  aday <- c(file.path(IN_DIR, "data", ad), file.path(IN_DIR, ad))
  var  <- aday[file.exists(aday)]
  if (!length(var)) stop(sprintf(
    "'%s' bulunamadi. Su konumlara bakildi:\n  %s\nDosyayi bunlardan birine koyun.",
    ad, paste(aday, collapse = "\n  ")))
  var[1]
}
e <- new.env(); load(bul("CALISMA2.Rdata"), envir = e)
d <- get(ls(e)[1], envir = e)

keys  <- readr::read_csv(bul("recovered_keys.csv"), show_col_types = FALSE)
items <- keys$madde_id
corr  <- setNames(keys$dogru_kod, keys$madde_id)
stopifnot(length(items) == length(unique(items)))
stopifnot(all(items %in% names(d)), all(paste0(items, "_R") %in% names(d)))

# Puanlama: yalniz 1-4 gecerli (6=ulasilamadi, 9=atlandi -> NA)
# Veri SPSS etiketli (haven_labelled) gelebilir; ham sayisal kodu al (etiketi soy).
zl <- function(x) {
  if (inherits(x, "haven_labelled") || inherits(x, "labelled")) x <- unclass(x)
  suppressWarnings(as.numeric(x))
}
id_all <- as.character(zl(d$IDSTUD))            # ogrenci kimligi (etiketsiz, karakter)
acc_df <- data.frame(IDSTUD = id_all, stringsAsFactors = FALSE)
for (it in items) { r <- zl(d[[it]])
  acc_df[[it]] <- ifelse(r %in% GECERLI, as.integer(r == corr[[it]]), NA_integer_) }

# Sure: ilk yanit suresi (_R), gecersizler -> NA (saniye)
fat_df <- data.frame(IDSTUD = id_all, stringsAsFactors = FALSE)
for (it in items) { t <- zl(d[[paste0(it, "_R")]])
  fat_df[[it]] <- ifelse(is.finite(t) & t > 0 & t < MAX_TIME, t, NA_real_) }

long <- acc_df |> tidyr::pivot_longer(-IDSTUD, names_to = "items", values_to = "acc") |>
  dplyr::left_join(
    tidyr::pivot_longer(fat_df, -IDSTUD, names_to = "items", values_to = "fat"),
    by = c("IDSTUD", "items"))
cat(sprintf("Ogrenci: %d | Madde: %d | Toplam gozlem: %d\n\n",
            length(unique(id_all)), length(items), nrow(long)))

# ===================== 2. NT10 + HIZLI TAHMIN SINIFLAMASI ===========
cat("== 2. NT10 (2-SS budamali) esigi ve hizli tahmin (RG) sinifi ==\n")
winsor_top <- function(x, k = WINSOR_SD) { ok <- is.finite(x) & x > 0
  m <- mean(x[ok]); s <- sd(x[ok]); cap <- m + k * s; ifelse(ok & x > cap, cap, x) }
nt_thr <- function(x, pct = .10, maxsec = 10) { x <- x[is.finite(x) & x > 0]
  if (length(x) < NT_MINN) return(NA_real_); min(pct * mean(x), maxsec) }

long$fatx <- ave(long$fat, long$items, FUN = function(x) winsor_top(x))
thr_tab <- long |> dplyr::filter(is.finite(fatx), fatx > 0) |>
  dplyr::group_by(items) |>
  dplyr::summarise(thr = nt_thr(fatx), .groups = "drop")
long <- long |> dplyr::left_join(thr_tab, by = "items") |>
  dplyr::mutate(rg = as.integer(!is.na(thr) & is.finite(fat) & fat > 0 & fat <= thr))

eng <- long |> dplyr::filter(is.finite(fat), fat > 0)
cat(sprintf("Genel RG orani: %.4f | caba-gosteren gozlem: %d\n\n",
            mean(eng$rg), nrow(eng)))

# --- RG gecerlik: hizli tahmin dogrulugu sansa (0.25) yakin olmali ---
val <- eng |> dplyr::filter(!is.na(acc)) |> dplyr::group_by(rg) |>
  dplyr::summarise(n = dplyr::n(), dogruluk = round(mean(acc), 3), .groups = "drop") |>
  dplyr::mutate(durum = ifelse(rg == 1, "Hizli tahmin", "Caba"))
cat("RG gecerlik kontrolu (sans = 0.25):\n"); print(as.data.frame(val)); cat("\n")
readr::write_csv(val, file.path(OUT_DIR, "T0_rg_gecerlik.csv"))

# ===================== 3. YANIT MATRISLERI =========================
# (a) Standart 2PL  : tum puanlanmis yanitlar (genis)
# (b) EM-IRT        : hizli tahmin yanitlari NA (Wise-DeMars ile ozdes kestirim)
# (c) IRTree        : pseudo-madde matrisi (item_n1 = RG, item_n2 = ACC | RG=0)
cat("== 3. Yanit matrisleri kuruluyor ==\n")
acc_wide <- acc_df; rownames(acc_wide) <- acc_wide$IDSTUD; acc_wide$IDSTUD <- NULL
acc_wide <- acc_wide[, items, drop = FALSE]

# EM-IRT: effortful yanitlarda acc, hizli tahminde NA
rg_wide <- long |> dplyr::select(IDSTUD, items, rg) |>
  tidyr::pivot_wider(names_from = items, values_from = rg)
rg_wide <- as.data.frame(rg_wide); rownames(rg_wide) <- rg_wide$IDSTUD
rg_wide <- rg_wide[rownames(acc_wide), items, drop = FALSE]
emirt_wide <- acc_wide
emirt_wide[which(rg_wide == 1, arr.ind = TRUE)] <- NA_integer_

# IRTree pseudo-madde matrisi
node <- eng |> dplyr::mutate(n1 = rg, n2 = ifelse(rg == 0, acc, NA_integer_))
n1w <- node |> dplyr::select(IDSTUD, items, n1) |>
  tidyr::pivot_wider(names_from = items, values_from = n1, names_glue = "{items}_n1")
n2w <- node |> dplyr::select(IDSTUD, items, n2) |>
  tidyr::pivot_wider(names_from = items, values_from = n2, names_glue = "{items}_n2")
it_in <- items[items %in% unique(node$items)]
irt_resp <- as.data.frame(dplyr::left_join(n1w, n2w, by = "IDSTUD"))
rownames(irt_resp) <- irt_resp$IDSTUD; irt_resp$IDSTUD <- NULL
irt_resp <- irt_resp[, c(paste0(it_in, "_n1"), paste0(it_in, "_n2")), drop = FALSE]
cat(sprintf("acc_wide: %d x %d | emirt_wide (NA'li) | irt_resp: %d x %d\n\n",
            nrow(acc_wide), ncol(acc_wide), nrow(irt_resp), ncol(irt_resp)))

# ===================== 4. MODELLER =================================
cat("== 4. Modeller kestiriliyor (surebilir) ==\n")

# (a) Standart 2PL -----------------------------------------------------
cat("  (a) Standart 2PL ...\n")
fit_std <- mirt::mirt(acc_wide, 1, itemtype = "2PL", verbose = FALSE)

# (b) EM-IRT (hizli tahmin -> NA, effortful'da 2PL) --------------------
cat("  (b) EM-IRT (Wise-DeMars) ...\n")
fit_em <- mirt::mirt(emirt_wide, 1, itemtype = "2PL", verbose = FALSE)

# (c) IRTree 2PL (RG + ACC dugum) -------------------------------------
cat("  (c) IRTree 2PL (RG + ACC) ...\n")
J <- length(it_in)
spec <- mirt::mirt.model(sprintf("RG = 1-%d\n ACC = %d-%d\n COV = RG*ACC",
                                  J, J + 1, 2 * J))
fit_tree <- mirt::mirt(irt_resp, spec, itemtype = "2PL", method = "EM",
                       technical = list(NCYCLES = 2000), verbose = FALSE)
cat("  Modeller tamam.\n\n")

# ===================== 5. MADDE PARAMETRELERI ======================
cat("== 5. Madde parametreleri toplaniyor (a, b) ==\n")
# Standart & EM-IRT: dogrudan IRTpars a,b
par_std <- mirt::coef(fit_std, simplify = TRUE, IRTpars = TRUE)$items
par_em  <- mirt::coef(fit_em,  simplify = TRUE, IRTpars = TRUE)$items
std_ab  <- data.frame(item = rownames(par_std), a_std = par_std[, "a"], b_std = par_std[, "b"])
em_ab   <- data.frame(item = rownames(par_em),  a_em  = par_em[, "a"],  b_em  = par_em[, "b"])

# IRTree: ACC dugumu (a = ACC boyutundaki yukleme; b = -d/a)
co_tree <- mirt::coef(fit_tree, simplify = TRUE)$items
acc_rows <- grepl("_n2$", rownames(co_tree))
tree_ab <- data.frame(
  item = sub("_n2$", "", rownames(co_tree)[acc_rows]),
  a_tree = co_tree[acc_rows, "a2"],
  b_tree = -co_tree[acc_rows, "d"] / co_tree[acc_rows, "a2"])

item_par <- std_ab |>
  dplyr::full_join(em_ab, by = "item") |>
  dplyr::full_join(tree_ab, by = "item")
readr::write_csv(item_par, file.path(OUT_DIR, "T1_madde_parametreleri.csv"))
cat("Madde parametreleri (ilk 6):\n"); print(head(item_par)); cat("\n")

# ===================== 6. YETENEK (THETA) ==========================
cat("== 6. Yetenek (theta) kestirimleri ==\n")
th_std <- mirt::fscores(fit_std, method = "EAP")[, 1]
th_em  <- mirt::fscores(fit_em,  method = "EAP")[, 1]
th_tree_mat <- mirt::fscores(fit_tree, method = "EAP")
th_tree <- th_tree_mat[, "ACC"]

theta <- data.frame(IDSTUD = rownames(acc_wide),
                    theta_std = th_std, theta_em = th_em)
theta_tree <- data.frame(IDSTUD = rownames(irt_resp), theta_tree = th_tree)
theta <- theta |> dplyr::left_join(theta_tree, by = "IDSTUD")
readr::write_csv(theta, file.path(OUT_DIR, "T2_yetenek_theta.csv"))

# ===================== 7. KARSILASTIRMA METRIKLERI =================
cat("== 7. Karsilastirma metrikleri ==\n")
karsilastir <- function(x, y, etiket) {
  ok <- is.finite(x) & is.finite(y)
  data.frame(karsilastirma = etiket, n = sum(ok),
    r_pearson = round(cor(x[ok], y[ok]), 3),
    rho_spearman = round(cor(x[ok], y[ok], method = "spearman"), 3),
    ort_fark = round(mean(y[ok] - x[ok]), 3),
    RMSD = round(sqrt(mean((y[ok] - x[ok])^2)), 3))
}
cmp_item <- dplyr::bind_rows(
  # gucluk b
  karsilastir(item_par$b_std, item_par$b_em,   "b: Standart -> EM-IRT"),
  karsilastir(item_par$b_std, item_par$b_tree, "b: Standart -> IRTree(ACC)"),
  karsilastir(item_par$b_em,  item_par$b_tree, "b: EM-IRT -> IRTree(ACC)"),
  # ayirt edicilik a
  karsilastir(item_par$a_std, item_par$a_em,   "a: Standart -> EM-IRT"),
  karsilastir(item_par$a_std, item_par$a_tree, "a: Standart -> IRTree(ACC)"),
  karsilastir(item_par$a_em,  item_par$a_tree, "a: EM-IRT -> IRTree(ACC)"))
cmp_theta <- dplyr::bind_rows(
  karsilastir(theta$theta_std, theta$theta_em,   "theta: Standart -> EM-IRT"),
  karsilastir(theta$theta_std, theta$theta_tree, "theta: Standart -> IRTree(ACC)"),
  karsilastir(theta$theta_em,  theta$theta_tree, "theta: EM-IRT -> IRTree(ACC)"))
cmp_all <- dplyr::bind_rows(cmp_item, cmp_theta)
readr::write_csv(cmp_all, file.path(OUT_DIR, "T3_karsilastirma_metrikleri.csv"))
cat("Karsilastirma metrikleri:\n"); print(cmp_all); cat("\n")

# ===================== 8. FIGURLER =================================
cat("== 8. Figurler ==\n")
sac <- function(df, xx, yy, xlab, ylab, baslik) {
  ggplot2::ggplot(df, ggplot2::aes(.data[[xx]], .data[[yy]])) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(color = "#1b7837", size = 2, alpha = .8) +
    ggplot2::labs(x = xlab, y = ylab, title = baslik) + theme_apa()
}
savefig(sac(item_par, "b_std", "b_em",   "Gucluk b - Standart", "Gucluk b - EM-IRT",
            "Sekil 1. Madde guclugu: Standart vs EM-IRT"), "fig1_b_std_em.png")
savefig(sac(item_par, "b_em", "b_tree",  "Gucluk b - EM-IRT", "Gucluk b - IRTree(ACC)",
            "Sekil 2. Madde guclugu: EM-IRT vs IRTree"), "fig2_b_em_tree.png")
savefig(sac(item_par, "a_std", "a_em",   "Ayirt edicilik a - Standart", "a - EM-IRT",
            "Sekil 3. Ayirt edicilik: Standart vs EM-IRT"), "fig3_a_std_em.png")
savefig(sac(theta,    "theta_em", "theta_tree", "theta - EM-IRT", "theta - IRTree(ACC)",
            "Sekil 4. Yetenek: EM-IRT vs IRTree"), "fig4_theta_em_tree.png")

# ===================== 8b. DUYARLILIK (esik kurallari) =============
# NT10 (budamali/budasiz), NT5, sabit 5 sn, MRTQ. Her kural icin RG orani,
# gecerlik ve ACC dugumu ort. ayirt ediciligi + ANA (NT10) ile korelasyon.
cat("== 8b. Duyarlilik analizi (esik kurallari) ==\n")

# MRTQ: log-RT 2-bilesenli normal karisim; hizli-tahmin bileseni agirligi w ->
# esik = exp( w-inci kuantil(log-RT) )  (Eker & Gelbal).
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

# Bir kurala gore uzun tabloya rg sutunu ekler
apply_rule <- function(L, method = "NT", pct = .10, maxsec = 10, fixed = NULL, winsor = TRUE) {
  L$fx <- if (winsor) ave(L$fat, L$items, FUN = function(x) winsor_top(x)) else L$fat
  th <- L |> dplyr::filter(is.finite(fx), fx > 0) |> dplyr::group_by(items) |>
    dplyr::summarise(thr = if (method == "MRTQ") mrtq_thr(fx)
                           else if (!is.null(fixed)) fixed
                           else nt_thr(fx, pct, maxsec), .groups = "drop")
  L |> dplyr::select(-dplyr::any_of("thr")) |> dplyr::left_join(th, by = "items") |>
    dplyr::mutate(rg = as.integer(!is.na(thr) & is.finite(fat) & fat > 0 & fat <= thr))
}

# Bir rg-sinifli uzun tablodan EM-IRT + IRTree kurup ozet doner
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
  eng2 <- Lr |> dplyr::filter(is.finite(fat), fat > 0, !is.na(acc))
  list(a_em = pe[, "a"], b_em = pe[, "b"], item = rownames(pe),
       theta_em = mirt::fscores(fe, method = "EAP")[, 1],
       rg_rate = mean(Lr$rg[is.finite(Lr$fat) & Lr$fat > 0]),
       acc_rg = mean(eng2$acc[eng2$rg == 1]), acc_sol = mean(eng2$acc[eng2$rg == 0]),
       a_acc_tree = mean(co[ar, "a2"]),
       rho = mirt::coef(ft, simplify = TRUE)$cov["ACC", "RG"] /
             sqrt(mirt::coef(ft, simplify = TRUE)$cov["ACC","ACC"] *
                  mirt::coef(ft, simplify = TRUE)$cov["RG","RG"]))
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
  cat("  -", nm, "...\n")
  args <- kurallar[[nm]]
  Lr <- do.call(apply_rule, c(list(L = Lbase), args))
  fr <- fit_rule(Lr)
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
cat("Duyarlilik tablosu:\n"); print(rob_tab); cat("\n")

# ===================== 9. KAYDET + OZET ============================
cat("== 9. Kaydediliyor ==\n")
saveRDS(list(item_par = item_par, theta = theta, cmp = cmp_all, val = val,
             thr = thr_tab, robustness = rob_tab,
             fits = list(std = fit_std, em = fit_em, tree = fit_tree)),
        file.path(OUT_DIR, "calisma3_TUM_SONUCLAR.rds"))

cat("\n##### OZET #####\n")
cat(sprintf("Ornek: %d ogrenci, %d madde | Genel RG: %.3f\n",
            nrow(acc_wide), length(items), mean(eng$rg)))
cat(sprintf("RG dogruluk: %.3f vs Caba dogruluk: %.3f (sans 0.25)\n",
            val$dogruluk[val$rg == 1], val$dogruluk[val$rg == 0]))
cat("Madde gucluk korelasyonlari (Standart/EM/IRTree) ve theta korelasyonlari icin T3'e bakin.\n")
cat("Cikti klasoru:", OUT_DIR, "\n")
cat("##### BITTI #####\n")
sink(); close(.logcon)
