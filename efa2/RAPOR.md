# Açımlayıcı Faktör Analizi Raporu — `ppass.xlsx`

**Veri:** 202 katılımcı, 21 madde, 5'li Likert
**Ters kodlanan maddeler:** M3, M4, M10, M11, M16

---

## 1. Sonuç: 3 faktörlü, 13 maddelik yapı

En iyi boyutluluk ve en iyi uyum iyiliği değerlerini **birlikte** veren yapı: **3 faktör, 13 madde.**

| Uyum indeksi | Değer | Ölçüt | Durum |
|---|---|---|---|
| χ² (robust) | 22,0 (sd = 22,6) | — | p > 0,05 |
| **RMSEA** | **0,000** | < 0,06 | ✔ |
| **CFI** | **1,000** | > 0,95 | ✔ |
| **TLI** | **1,004** | > 0,95 | ✔ |
| **SRMR** | **0,056** | < 0,08 | ✔ |
| KMO | 0,743 | > 0,70 | ✔ |
| En düşük MSA | 0,704 | > 0,50 | ✔ |
| Bartlett | χ² = 408,5; sd = 78; p < 0,001 | p < 0,05 | ✔ |
| Açıklanan varyans | %42,3 | — | — |
| Faktörler arası r | 0,31 – 0,32 | < 0,85 | ✔ ayırt edici |
| Binişik yük | **0** | — | ✔ |
| En düşük faktör yükü | **0,44** | > 0,40 | ✔ |

### Faktör deseni

| Faktör | Madde | Yük | F1 | F2 | F3 | h² | MSA | Ort. | SS |
|---|---|---|---|---|---|---|---|---|---|
| **F1** | M7 | **0,650** | 0,65 | 0,01 | −0,02 | 0,417 | 0,704 | 3,80 | 1,37 |
| | M8 | **0,584** | 0,58 | 0,17 | 0,01 | 0,433 | 0,745 | 3,95 | 1,31 |
| | M1 | **0,568** | 0,57 | 0,09 | −0,07 | 0,337 | 0,706 | 3,63 | 1,34 |
| | M2 | **0,536** | 0,54 | −0,13 | 0,19 | 0,347 | 0,706 | 4,36 | 1,06 |
| | M21 | **0,489** | 0,49 | 0,00 | 0,29 | 0,415 | 0,813 | 4,43 | 0,90 |
| **F2** | M3ᴿ | **0,729** | 0,19 | 0,73 | 0,03 | 0,669 | 0,766 | 4,56 | 1,02 |
| | M10ᴿ | **0,697** | −0,03 | 0,70 | 0,05 | 0,495 | 0,727 | 4,72 | 0,86 |
| | M4ᴿ | **0,622** | 0,33 | 0,62 | −0,13 | 0,566 | 0,779 | 4,60 | 1,03 |
| | M11ᴿ | **0,450** | −0,07 | 0,45 | 0,17 | 0,255 | 0,765 | 4,35 | 1,19 |
| **F3** | M19 | **0,669** | 0,09 | 0,06 | 0,67 | 0,525 | 0,761 | 4,37 | 1,03 |
| | M20 | **0,666** | −0,13 | 0,27 | 0,67 | 0,570 | 0,735 | 4,12 | 1,13 |
| | M9 | **0,453** | 0,06 | −0,10 | 0,45 | 0,206 | 0,716 | 4,72 | 0,67 |
| | M18 | **0,439** | 0,06 | 0,13 | 0,44 | 0,266 | 0,799 | 3,30 | 1,21 |

ᴿ = ters kodlanmış madde

### Güvenirlik

| Faktör | Madde | Sıralı α | McDonald ω | Ort. yük |
|---|---|---|---|---|
| F1 | M1, M2, M7, M8, M21 | 0,741 | 0,724 | 0,57 |
| F2 | M3, M4, M10, M11 | 0,763 | 0,756 | 0,62 |
| F3 | M9, M18, M19, M20 | 0,687 | 0,671 | 0,56 |
| **Toplam** | 13 madde | **0,813** | — | — |

### Çıkarılan 8 madde
`M5, M6, M12, M13, M14, M15, M16, M17` — düşük ortak varyans (h² < 0,30) veya
düşük faktör yükü (< 0,40) nedeniyle, teker teker ve her adımda çözüm yeniden
kestirilerek çıkarılmıştır. En sorunlu olanlar: **M12** (h² = 0,09–0,15;
madde-toplam r = 0,18) ve **M17** (MSA = 0,58; h² = 0,08–0,32).

---

## 2. Faktör Sayısı Neden 3?

| Yöntem | Sonuç |
|---|---|
| **Paralel analiz — PC, ortalama** (500 permütasyon) | **3** |
| **Paralel analiz — PC, %95** | **3** |
| **Paralel analiz — indirgenmiş matris, %95** | **3** |
| Paralel analiz — indirgenmiş matris, ortalama | 4 |
| Velicer MAP (4. kuvvet) | 3 |
| Velicer MAP (kare) | 1 |
| Kaiser (özdeğer > 1) | 7 |
| Ampirik Kaiser Ölçütü | 7 |
| 1./2. özdeğer oranı | 2,87 (< 4 → **tek boyutlu değil**) |

Kaiser ölçütü ve Ampirik Kaiser Ölçütü bilindiği gibi aşırı faktörleştirmektedir
(21 maddede 7 faktör). En güvenilir yöntem olan paralel analiz dört ölçütün üçünde
**3** vermektedir.

### Uyum karşılaştırması (robust WLSMV tipi)

| Madde havuzu | k | RMSEA | CFI | TLI | SRMR | En düşük yük | h²<0,30 | Binişik | Φ_maks |
|---|---|---|---|---|---|---|---|---|---|
| 21 (tümü) | 1 | 0,055 | 0,839 | 0,878 | 0,105 | 0,23 | 16 | 0 | — |
| 21 (tümü) | 2 | 0,042 | 0,915 | 0,930 | 0,087 | 0,28 | 11 | 0 | 0,39 |
| 21 (tümü) | 3 | 0,025 | 0,967 | 0,974 | 0,069 | 0,32 | 9 | 1 | 0,33 |
| 21 (tümü) | 4 | 0,000 | 1,000 | 1,006 | 0,057 | 0,28 | 4 | 2 | 0,31 |
| 16 madde | 3 | 0,000 | 1,000 | 1,008 | 0,060 | 0,37 | 7 | 0 | 0,32 |
| 15 madde | 3 | 0,000 | 1,000 | 1,015 | 0,056 | 0,32 | 5 | 0 | 0,36 |
| **13 madde** | **3** | **0,000** | **1,000** | **1,004** | **0,056** | **0,44** | **3** | **0** | **0,32** |
| 13 madde | 1 | 0,067 | 0,821 | 0,849 | 0,113 | 0,28 | 9 | 0 | — |
| 13 madde | 2 | 0,055 | 0,890 | 0,900 | 0,088 | 0,17 | 4 | 5 | 0,47 |

Birden fazla model mükemmel uyum vermektedir; 13 maddelik 3 faktörlü çözüm bunlar
arasında **en yüksek en-düşük-yüke (0,44), en az düşük ortak varyanslı maddeye ve
sıfır binişik yüke** sahip olanıdır.

### Yapı kararlılığı — ve neden ham değeri okumak yanıltıcı olurdu

Yarı-örneklem çapraz geçerleme (150 tekrar, Tucker uyum katsayısı φ) ham hâlde
k = 1'i işaret ediyor gibi görünmektedir. Ancak **n = 202'de yarılar yalnızca
n ≈ 101'dir** ve bu, küçük faktörler için son derece sert bir testtir. Bu nedenle
gözlenen yükler ve eşikler kullanılarak, **gerçekten 3 faktörlü** bir modelden
n = 202 ile veri üretilip aynı prosedür uygulanmış ve bir referans tavanı
hesaplanmıştır:

| k | Gözlenen φ | Referans φ (yapı gerçekten 3 faktörlüyken) | Değerlendirme |
|---|---|---|---|
| 1 | 0,934 | 0,929 | tavanda |
| 2 | 0,592 | 0,740 | **tavanın belirgin altında → kararsız** |
| **3** | **0,696** | **0,731** | **tavanda (%95'i) → beklendiği kadar kararlı** |
| 4 | 0,595 | 0,648 | tavanın altında |

Yani bu örneklem büyüklüğünde φ ≈ 0,73, doğru bir 3 faktörlü yapı için ulaşılabilecek
**en yüksek** değerdir; gözlenen 0,696 pratikte bu tavandadır. Buna karşılık 2 faktörlü
çözüm kendi referansının açıkça altında kalmaktadır. Kalibrasyon yapılmasaydı,
3 faktörlü yapı yanlışlıkla "kararsız" ilan edilecekti.

---

## 3. Yöntem

- **Polikorik korelasyon matrisi.** 5'li sıralı ve çarpık veride Pearson ilişkileri
  sistematik olarak zayıflatır. Tahminci bilinen doğruya karşı doğrulanmıştır (hata ±0,03).
  Seyrek hücreli 4 maddede (M4, M9, M10, M21) bitişik kategoriler birleştirilmiştir.
- **Çıkarım:** MINRES/ULS. **Döndürme:** oblimin (eğik), Kaiser normalizasyonlu,
  Gradient Projection Algorithm (Jennrich, 2002) ile.
- **Robust uyum indeksleri.** Polikorik matris üzerinde naif ML ki-karesi ciddi
  biçimde şişkindir — bu veride naif RMSEA 0,110 iken robust karşılığı 0,000'dır.
  Simülasyonla doğrulanmıştır: gerçekten 3 faktörlü, ağır çarpık veride naif ML doğru
  modele RMSEA = 0,098 / CFI = 0,893 ("kötü uyum") derken, robust düzeltme
  RMSEA = 0,029 / CFI = 0,975 vermektedir. Raporlanan indeksler, polikorik
  korelasyonların asimptotik kovaryans matrisinden (Γ; 1000 bootstrap) türetilen
  ortalama-ve-varyans düzeltmeli (WLSMV tipi) ki-kareye dayanır
  (Muthén, 1984; Satorra & Bentler, 1994; Asparouhov & Muthén, 2010).
- **Sistematik arama.** k ∈ {2,3,4,5} × h² eşiği ∈ {yok, 0,20, 0,25, 0,30} ×
  yük eşiği ∈ {0,32, 0,40, 0,45} × binişiklik farkı ∈ {0,15, 0,20} ızgarasında
  toplam 96 arındırma yolu denenmiş; her faktöründe en az 3 madde bulunan 11 geçerli
  çözüm elde edilmiş ve bunlar uyum, yük kalitesi, ayırt edicilik ve kararlılık
  açısından karşılaştırılmıştır. k = 4 ve k = 5 için geçerli çözüm üretilememiştir
  (her zaman 3 maddeden az faktör oluşmuştur).

> **Teknik not.** `factor_analyzer` kütüphanesinin eğik döndürmeleri k ≥ 3 için desen
> matrisiyle tutarsız bir faktör korelasyon matrisi döndürmekte
> (ΛΛ′ = ΛΦΛ′ değişmezliği 0,73'e varan farklarla ihlal edilmekte) ve uyum
> hesaplarını bozmaktadır. Döndürme bu nedenle bağımsız olarak uygulanmış,
> değişmezlik makine hassasiyetinde (4×10⁻¹⁶) doğrulanmıştır.

### Veri tarama
- Metin biçiminde kaydedilmiş 4 sütunda (M4, M6, M15, M16) `' '` eksik-veri işareti
  sayısala çevrilmiştir. Toplam eksik: 34 hücre (%0,80); ikili tam gözlem kullanılmıştır.
- Ters kodlama sonrası **hiçbir maddede negatif madde-toplam korelasyonu kalmamıştır**
  (M3: −0,219 → +0,407; M4: −0,185 → +0,364; M10: −0,047 → +0,227;
  M11: −0,074 → +0,199; M16: −0,224 → +0,276), yani belirtilen 5 maddenin ters
  kodlanması gerektiği veriyle doğrulanmıştır.

---

## 4. Önemli Uyarılar

**1. F2 bir ifade-yönü (method) faktörü olabilir.** F2'yi oluşturan dört maddenin
(**M3, M4, M10, M11**) tamamı ters kodlanmış maddedir. Ters ifade edilmiş maddelerin
kendi aralarında bir faktör oluşturması, ölçme literatüründe iyi bilinen bir
**yöntem artefaktıdır** ve içeriksel bir alt boyutla karıştırılmamalıdır. Bu ayrımı
yapabilmek için madde metinleri gereklidir:

- F2 maddelerinin içeriği birbirine benziyorsa → gerçek bir alt boyut.
- İçerikleri birbirinden farklı, ortak yanları yalnızca olumsuz ifade edilmiş olmaları ise
  → yöntem faktörü; bu durumda yapı 2 faktörlü olarak yeniden değerlendirilmelidir.

Madde metinlerini paylaşırsanız bunu doğrudan sınayabilirim.

**2. Güvenirlikler orta düzeydedir.** F3 için α = 0,687, 0,70 eşiğinin hemen altındadır.
Faktörler 4–5 maddeli olduğundan bu beklenen bir durumdur, ancak F3 için madde
eklenmesi yararlı olur.

**3. Zayıf maddeler.** M9 (h² = 0,206) ve M11 (h² = 0,255) 0,30 eşiğinin altındadır;
faktörlerini 3 maddeye düşürmemek için korunmuşlardır.

**4. Ölçek genel olarak zayıftır.** Birinci faktör varyansın yalnızca %26,3'ünü
açıklamakta, 21 maddenin 7'sinde madde-toplam korelasyonu 0,30'un altında kalmaktadır.
21 maddenin 8'inin çıkarılması bunun sonucudur.

**5. Örneklem.** n = 202, 13 madde için madde başına ~16 katılımcıdır (yeterli), ancak
21 madde için sınırdadır. Bulgular **açımlayıcıdır**; yapının **bağımsız bir örneklemde
DFA** ile sınanması gerekir.

---

## 5. Alternatif Modeller

Daha fazla madde korunması isteniyorsa, uyum yine mükemmeldir:

| Model | Madde | k | RMSEA | CFI | SRMR | α | Ödün |
|---|---|---|---|---|---|---|---|
| **Önerilen** | **13** | **3** | **0,000** | **1,000** | **0,056** | **0,813** | 8 madde kaybı |
| Alternatif B | 15 | 3 | 0,000 | 1,000 | 0,056 | 0,830 | en düşük yük 0,32; 5 madde h² < 0,30 |
| Alternatif A | 16 | 3 | 0,000 | 1,000 | 0,060 | 0,815 | en düşük yük 0,37; 7 madde h² < 0,30 |
| Madde atmadan | 21 | 3 | 0,025 | 0,967 | 0,069 | 0,851 | 9 madde h² < 0,30; 1 binişik |

Üçünde de faktör bileşimi büyük ölçüde aynıdır; fark, zayıf maddelerin korunup
korunmamasındadır.

---

## 6. Dosyalar

| Dosya | İçerik |
|---|---|
| `01_prep.py` | Temizleme, ters kodlama, veri tarama |
| `02_gamma.py` | Bootstrap Γ (1000 tekrar) |
| `03_retain.py` | Özdeğer, MAP, EKC, paralel analiz |
| `05_search.py`, `06_grid.py` | Sistematik arındırma ızgarası |
| `08_stab.py` | Yarı-örneklem kararlılık |
| `10_calib.py` | Kararlılık referans tavanı (simülasyon) |
| `09_robust.py`, `11_final.py` | Robust uyum ve final çözüm |
| `final_13madde_3faktor.csv` | Final madde istatistikleri |
| `polychor.py`, `rotate.py`, `efa_core.py`, `robust_fit.py` | Yöntem altyapısı |

## Kaynaklar

Asparouhov, T., & Muthén, B. (2010). *Simple second order chi-square correction*. Mplus Technical Appendix.
Braeken, J., & van Assen, M. A. L. M. (2017). An empirical Kaiser criterion. *Psychological Methods, 22*(3), 450–466.
Flora, D. B., & Curran, P. J. (2004). An empirical evaluation of alternative methods of estimation for CFA with ordinal data. *Psychological Methods, 9*(4), 466–491.
Horn, J. L. (1965). A rationale and test for the number of factors in factor analysis. *Psychometrika, 30*(2), 179–185.
Hu, L., & Bentler, P. M. (1999). Cutoff criteria for fit indexes in covariance structure analysis. *Structural Equation Modeling, 6*(1), 1–55.
Jennrich, R. I. (2002). A simple general method for oblique rotation. *Psychometrika, 67*(1), 7–19.
Lorenzo-Seva, U., & ten Berge, J. M. F. (2006). Tucker's congruence coefficient as a meaningful index of factor similarity. *Methodology, 2*(2), 57–64.
Muthén, B. (1984). A general structural equation model with dichotomous, ordered categorical, and continuous latent variable indicators. *Psychometrika, 49*(1), 115–132.
Olsson, U. (1979). Maximum likelihood estimation of the polychoric correlation coefficient. *Psychometrika, 44*(4), 443–460.
Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics and standard errors in covariance structure analysis. In *Latent variables analysis* (pp. 399–419). Sage.
Velicer, W. F. (1976). Determining the number of components from the matrix of partial correlations. *Psychometrika, 41*(3), 321–327.
