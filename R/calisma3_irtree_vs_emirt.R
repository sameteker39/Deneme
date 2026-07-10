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
gerekli <- c("dplyr", "tidyr", "ggplot2", "mirt", "readr")
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
NA_TIME   <- 1e10          # gecersiz sure kodu
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
e <- new.env(); load(file.path(IN_DIR, "data", "CALISMA2.Rdata"), envir = e)
d <- get(ls(e)[1], envir = e)

keys  <- readr::read_csv(file.path(IN_DIR, "data", "recovered_keys.csv"),
                         show_col_types = FALSE)
items <- keys$madde_id
corr  <- setNames(keys$dogru_kod, keys$madde_id)
stopifnot(length(items) == length(unique(items)))
stopifnot(all(items %in% names(d)), all(paste0(items, "_R") %in% names(d)))

# Puanlama: yalniz 1-4 gecerli (6=ulasilamadi, 9=atlandi -> NA)
zl <- function(x) suppressWarnings(as.numeric(x))
acc_df <- data.frame(IDSTUD = d$IDSTUD)
for (it in items) { r <- zl(d[[it]])
  acc_df[[it]] <- ifelse(r %in% GECERLI, as.integer(r == corr[[it]]), NA_integer_) }

# Sure: ilk yanit suresi (_R), gecersizler -> NA (saniye)
fat_df <- data.frame(IDSTUD = d$IDSTUD)
for (it in items) { t <- zl(d[[paste0(it, "_R")]])
  fat_df[[it]] <- ifelse(is.finite(t) & t > 0 & t < NA_TIME, t, NA_real_) }

long <- acc_df |> tidyr::pivot_longer(-IDSTUD, names_to = "items", values_to = "acc") |>
  dplyr::left_join(
    tidyr::pivot_longer(fat_df, -IDSTUD, names_to = "items", values_to = "fat"),
    by = c("IDSTUD", "items"))
cat(sprintf("Ogrenci: %d | Madde: %d | Toplam gozlem: %d\n\n",
            length(unique(d$IDSTUD)), length(items), nrow(long)))

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

# ===================== 9. KAYDET + OZET ============================
cat("== 9. Kaydediliyor ==\n")
saveRDS(list(item_par = item_par, theta = theta, cmp = cmp_all, val = val,
             thr = thr_tab, fits = list(std = fit_std, em = fit_em, tree = fit_tree)),
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
