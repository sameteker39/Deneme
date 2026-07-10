# ÇALIŞMA 3 — IRTree ve EM-IRT Karşılaştırması (TIMSS 2023 Türkiye)

Bu depo, **hızlı tahmin (rapid guessing)** dikkate alındığında **madde ve yetenek
parametrelerinin** iki farklı yaklaşımla nasıl kestirildiğini karşılaştırır:

- **IRTree** (Item Response Tree) — iki düğümlü ağaç: `n1` = hızlı tahmin (RG),
  `n2` = doğruluk (ACC); doğruluk düğümü, RG boyutuyla **eş zamanlı** kestirilir.
- **EM-IRT** (Effort-Moderated IRT; Wise & DeMars, 2006) — hızlı tahmin yanıtları
  şans fonksiyonuna sabitlenir, kalan (effortful) yanıtlarda IRT.
- **Standart 2PL** — çabayı yok sayan temel model (baseline).

## Veri

- **Kaynak:** TIMSS 2023, Türkiye (`IDCNTRY = 792`), 8. sınıf matematik.
- **Örneklem:** 4925 öğrenci.
- **Maddeler:** 43 çoktan seçmeli (4 seçenekli) matematik maddesi. Her madde için
  hem ham yanıt (`ME#####`, 1–4) hem **ilk yanıt süresi** (`ME#####_R`, saniye) mevcut.
- **Ham veri dosyası** (`data/CALISMA2.Rdata`) bu depoya **dahil edilmemiştir**
  (bkz. `.gitignore`); analiz için `data/` altına konmalıdır.

### Cevap anahtarı (`data/recovered_keys.csv`)

43 maddenin cevap anahtarı, verideki TIMSS başarı puanları (plausible values,
`BSMMAT01–05`) kullanılarak **ampirik** olarak geri kazanılmıştır: her maddede doğru
seçenek, onu işaretleyen öğrencilerin ortalama başarısı en yüksek olan seçenektir.

- Gerçek anahtarı bilinen **5 maddede 5/5 doğrulanmıştır** (klasik distraktör analizi).
- 42 maddede seçenekler arası başarı farkı çok büyüktür (yüksek güven).
- Tek istisna **`ME72022`** (fark ≈ 15.7 puan) — resmi TIMSS anahtarıyla teyit edilmelidir.

> Elinizde resmi `madde_bilgi.xlsx` (`dogru_kod`) varsa, `recovered_keys.csv` yerine
> onu kullanabilir ya da çapraz kontrol edebilirsiniz.

## Yöntem

- **Eşik (rapid guessing):** NT10 (Wise & Ma, 2012) — madde başına
  `eşik = min(0.10 × ortalama_süre, 10 sn)`. NT hesabından önce uç süreler
  madde içinde **ortalama + 2 SS**'te budanır (winsorize; Goldhammer vd., 2014;
  Nagy, Ulitzsch & Lindner, 2019). Süresi eşiğin altındaki yanıt = hızlı tahmin.
- **Modeller** `mirt` ile kestirilir; hepsi 2PL, ortak gizil metrik.
- **Karşılaştırma:** madde güçlük (`b`), ayırt edicilik (`a`) ve yetenek (`theta`)
  için Pearson/Spearman korelasyonu, ortalama fark ve RMSD; saçılım grafikleri.

### Önemli yöntem notu — EM-IRT ↔ effortful-only özdeşliği

Wise & DeMars EM-IRT modelinde hızlı tahmin yanıtları şans fonksiyonuna sabitlenir;
bu yanıtların olabilirliğe katkısı **parametreden bağımsız bir sabittir**. Dolayısıyla
madde ve yetenek parametrelerinin kestirimi açısından bu model, "hızlı tahmin
yanıtlarını **NA** yapıp kalan yanıtlarda 2PL kestirmek" ile **özdeştir**. Betik bu
özdeş kestirimi kullanır. Fark, test bilgisi / skorlama yorumunda ortaya çıkar.

Bu nedenle **EM-IRT** ile **IRTree ACC düğümü** aynı effortful yanıt kümesini kullanır
(bu veride 26.355 gözlem). Aralarındaki fark, IRTree'nin RG boyutuyla **çok boyutlu
eş zamanlı** kestiriminden kaynaklanır — çalışmanın odak noktası budur.

## Çalıştırma

```r
# R (mirt, dplyr, tidyr, ggplot2, readr kurulu olmalı)
setwd("bu/deponun/yolu")
# data/CALISMA2.Rdata ve data/recovered_keys.csv yerinde olmalı
source("R/calisma3_irtree_vs_emirt.R")
```

Çıktılar `outputs/` altına yazılır:

| Dosya | İçerik |
|---|---|
| `T0_rg_gecerlik.csv` | RG geçerlik (doğruluk × durum) |
| `T1_madde_parametreleri.csv` | `a, b` (Standart / EM-IRT / IRTree-ACC) |
| `T2_yetenek_theta.csv` | Öğrenci `theta` kestirimleri (3 model) |
| `T3_karsilastirma_metrikleri.csv` | Korelasyon, ort. fark, RMSD |
| `fig1–4_*.png` | Saçılım grafikleri (600 dpi) |
| `calisma3_TUM_SONUCLAR.rds` | Tüm nesneler + model uyumları |

## Ön bulgular (doğrulama koşusundan)

Veri hazırlama + NT10 + sınıflama adımları gerçek veride doğrulanmıştır:

- Genel hızlı tahmin oranı ≈ **%11.2** (30.212 çaba-gösteren gözlem).
- RG geçerliği: hızlı tahmin doğruluğu ≈ **0.399**, çaba doğruluğu ≈ **0.550**
  (şans = 0.25) — beklenen yönde ayrışma.
- Öğrenci başı medyan 6 madde (kitapçık/matris tasarımı): madde parametreleri
  sağlam (madde başı ~700 kişi); yetenek `theta` daha az kesin (sınırlılık).

## Ön sonuçlar (Python, `mirt` olmadan)

`mirt` kurulamayan ortamlar için `python/` altında bağımsız bir MML-EM kestirici
vardır (`irt_est.py`; 2PL ve iki-boyutlu IRTree, eksik veri destekli, simülasyonla
doğrulanmış). `python/run_preliminary.py` üç modeli çalıştırıp `outputs_python/`
altına karşılaştırma tablolarını yazar. **Bu sonuçlar ön/onay amaçlıdır; nihai
analiz `mirt` betiğiyle yapılmalıdır.**

> **Süre sentinel'i:** Eksik yanıt süresi verideki `_R` sütunlarında `9999999999`
> ile kodlanmıştır. `t < 1e10` filtresi bunu **dışlamaz** (9999999999 < 1e10);
> geçerli üst sınır `1e5` alınmıştır (gerçek RT < ~900 sn). Bu düzeltme yapılmazsa
> NT ortalaması şişer ve eşik 10 sn'de tıkanır.

Bu veride bulunan başlıca ön bulgular (NT10, RG oranı ≈ %6.2):

- Ortalama ayırt edicilik: Standart **1.430** → EM-IRT **1.449** → IRTree-ACC **1.456**
  (çabayı yok saymak ayırt ediciliği hafifçe zayıflatıyor).
- **EM-IRT ≈ IRTree-ACC**: madde güçlüğü `b` için r = **1.00**, ayırt edicilik `a`
  için r = **0.999** (özdeşlik ampirik olarak doğrulandı).
- Standart → EM-IRT: `b` için RMSD = **0.06**, `a` için RMSD = **0.10** (ılımlı ama
  tutarlı yönde); en çok kayan zor madde **ME62341** (b: 1.80 → 2.07).
- **IRTree gizil korelasyonu RG × ACC ≈ −0.32**: hızlı tahmin eğilimi düşük
  matematik yeteneğiyle ilişkili (IRTree'nin EM-IRT üzerine kişi düzeyindeki katkısı).

### Duyarlılık analizi (`python/run_robustness.py` → `outputs_python/T5_duyarlilik.csv`)

Beş eşik kuralı: NT10 (budamalı/budasız), NT5, sabit 5 sn, **MRTQ** (log-RT 2-bileşenli
karışım kuantili; Eker & Gelbal). RG oranı %3.3 (NT5) – %14.2 (MRTQ) aralığında değişir.
**Madde ve yetenek parametreleri eşik seçimine karşı çok sağlamdır**: EM-IRT `a`/`b`/`θ`
kestirimlerinin ANA (NT10) ile korelasyonu tüm kurallarda **≥ 0.99**. Eşiğe duyarlı tek
parametre, RG–ACC gizil korelasyonudur (NT5'te −0.15, MRTQ'da −0.42) — beklenen biçimde,
hızlı tahminin ne kadar yakalandığına bağlıdır.

## Kaynaklar

- Wise, S. L., & DeMars, C. E. (2006). *An application of item response time: The
  effort-moderated IRT model.* Journal of Educational Measurement, 43(1), 19–38.
- Wise, S. L., & Ma, L. (2012). *Setting response time thresholds for a CAT item
  pool: The normative threshold method.* (NCME).
- De Boeck, P., & Partchev, I. (2012). *IRTrees.* Journal of Statistical Software, 48.
