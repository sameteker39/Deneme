# Açımlayıcı Faktör Analizi Raporu

**Veri:** `FAfA_veri_31madde.xlsx` — 310 katılımcı, 31 madde, 5'li Likert
**Analiz tarihi:** 2026-08-06

---

## 1. Yönetici Özeti

31 maddelik havuza uygulanan açımlayıcı faktör analizi sonucunda, **9 madde literatüre uygun
ölçütlerle çıkarılmış** ve geriye kalan **22 maddenin tek boyutlu (unidimensional) bir yapı**
oluşturduğu belirlenmiştir.

| | Değer | Ölçüt | Durum |
|---|---|---|---|
| Madde sayısı | 22 (31'den) | — | 9 madde atıldı |
| Faktör sayısı | **1** | — | — |
| Açıklanan varyans | %45,4 | — | — |
| RMSEA | **0,048** | < 0,06 iyi | ✔ |
| CFI | **0,969** | > 0,95 iyi | ✔ |
| TLI | **0,982** | > 0,95 iyi | ✔ |
| SRMR | **0,077** | < 0,08 kabul | ✔ |
| KMO | 0,899 | > 0,80 çok iyi | ✔ |
| Sıralı α / ω | 0,947 / 0,948 | > 0,70 | ✔ |
| Faktör yükleri | 0,52 – 0,82 | > 0,40 | ✔ |
| **Yapı kararlılığı (Tucker φ)** | **0,993** | > 0,95 eşdeğer | ✔ |

**Kritik bulgu:** Çok faktörlü çözümler (k = 2, 3, 4) sayısal olarak *daha iyi* uyum indeksleri
üretmektedir (ör. k = 3 için RMSEA = 0,023, CFI = 0,994). Ancak yarı-örneklem çapraz geçerleme,
bu faktörlerin **replike olmadığını** göstermiştir (Tucker φ = 0,65–0,75; en zayıf faktör 0,54).
Yani çok faktörlü çözümlerin uyum üstünlüğü **şans eseri varyansın kapitalize edilmesinden**
kaynaklanmaktadır. Tek faktörlü çözüm ise neredeyse kusursuz replike olmaktadır (φ = 0,993).

---

## 2. Ön Analizler ve Veri Tarama

### 2.1 Eksik veri
64 hücre eksik (%0,67), 49 satırda dağılmış durumda. Oran çok düşük olduğu için liste bazlı
silme (n = 261'e düşerdi) yerine **ikili tam gözlem (pairwise)** yaklaşımı kullanılmıştır
(ikili n: 296–310, ortalama 306).

### 2.2 Ters kodlu madde — M9
M9, diğer 30 maddenin 27'siyle **negatif** ilişkiliydi (ortalama r = −0,22; düzeltilmiş
madde-toplam korelasyonu = **−0,47**). Bu, ters ifade edilmiş ancak ters kodlanmamış bir maddenin
tipik imzasıdır. M9 yeniden kodlanmıştır (6 − x); sonrasında madde-toplam korelasyonu **+0,47**
olmuştur.

> ⚠️ Bu düzeltme yapılmadan yürütülen herhangi bir analiz hatalı sonuç verir. Veri setinde M9'un
> ham (ters kodlanmamış) hâlde bulunduğunu doğrulamanız önerilir.

### 2.3 Normallik
Dağılımlar ciddi biçimde çarpıktır: 15 maddede |çarpıklık| > 2, bazı maddelerde basıklık 20'ye
ulaşmaktadır (M18: −4,44 / 20,32). Bu, Likert ölçeklerinde yaygın olan güçlü **tavan etkisidir**
ve yöntem seçimini doğrudan belirlemiştir (bkz. Bölüm 3).

### 2.4 Örneklem uygunluğu
| | 31 madde | 22 madde (final) |
|---|---|---|
| KMO | 0,852 | **0,899** |
| Bartlett χ² | 3529,98 (p < 0,001) | 2360,42 (p < 0,001) |
| En düşük madde MSA | 0,567 (M11) | **0,664** |

Her iki aşamada da veri faktör analizine uygundur; arındırma sonrası belirgin iyileşme vardır.

---

## 3. Yöntemsel Kararlar ve Gerekçeleri

### 3.1 Polikorik korelasyon kullanımı
5'li sıralı ve ağır çarpık veride Pearson korelasyonu ilişkileri sistematik olarak **zayıflatır**.
Simülasyonla doğrulanmıştır: gerçek r = 0,50 olan aşırı çarpık iki madde için Pearson 0,25,
polikorik 0,45 vermektedir. Gerçek veride polikorik ile Pearson arasındaki ortalama mutlak fark
**0,149**'dur — göz ardı edilemeyecek bir büyüklük. Bu nedenle tüm analizler **polikorik korelasyon
matrisi** üzerinde yürütülmüştür (Flora & Curran, 2004; Holgado-Tello vd., 2010).

Polikorik tahminci, bilinen doğruya karşı doğrulanmıştır (hata ±0,03; iki değişkenli normal
dağılım fonksiyonu makine hassasiyetinde).

### 3.2 Robust (WLSMV tipi) uyum indeksleri — kritik düzeltme
Polikorik matris üzerinde **naif ML ki-karesi hesaplamak ciddi biçimde yanıltıcıdır.**
Bu, doğrudan simülasyonla gösterilmiştir:

> Gerçekten 3 faktörlü, n = 310, ağır tavan etkili yapay veri üretildi.
> **Doğru modelde (k = 3):** naif ML → RMSEA = 0,098 / CFI = 0,893 (yani "kötü uyum" yargısı,
> model tam doğru olmasına rağmen). Robust düzeltme → RMSEA = 0,029 / CFI = 0,975 (doğru yargı).

Bu nedenle raporlanan uyum indeksleri, polikorik korelasyonların **asimptotik kovaryans
matrisinden (Γ, 1000 bootstrap ile kestirildi)** türetilen ortalama-ve-varyans düzeltmeli
(WLSMV tipi) ki-kare üzerine kuruludur (Muthén, 1984; Satorra & Bentler, 1994;
Asparouhov & Muthén, 2010).

Karşılaştırma için final modelde: naif ML → RMSEA = 0,154 / CFI = 0,676; robust → RMSEA = 0,048 /
CFI = 0,969.

### 3.3 Çıkarım ve döndürme
- **Çıkarım:** MINRES / ULS (ağır normallik ihlali nedeniyle ML yerine)
- **Döndürme:** Oblimin (eğik), Kaiser normalizasyonlu, Gradient Projection Algorithm
  (Jennrich, 2002) ile

> **Teknik not:** Yaygın kullanılan `factor_analyzer` kütüphanesinin eğik döndürmeleri
> (oblimin/promax/quartimin), k ≥ 3 için desen matrisiyle **tutarsız** bir faktör korelasyon
> matrisi (Φ) döndürmektedir — ΛΛ′ = ΛΦΛ′ değişmezliği 0,73'e varan farklarla ihlal edilmektedir.
> Bu, uyum indekslerini tamamen bozar. Analizde GPA döndürmesi bağımsız olarak uygulanmış ve
> değişmezlik makine hassasiyetinde (4×10⁻¹⁶) doğrulanmıştır.

---

## 4. Faktör Sayısının Belirlenmesi

### 4.1 Başlangıç (31 madde)

| Yöntem | Önerilen k |
|---|---|
| Kaiser (özdeğer > 1) | 7 |
| Ampirik Kaiser Ölçütü (Braeken & van Assen, 2017) | 7 |
| Velicer MAP (kare) | 4 |
| Velicer MAP (4. kuvvet) | 5 |
| **Paralel analiz — PC (500 permütasyon)** | **3** |
| **Paralel analiz — indirgenmiş matris, %95** | **4** |

Kaiser ölçütü bilindiği üzere aşırı faktörleştirmektedir; yakınsayan kanıt 3–5 faktöre
işaret etmiştir.

### 4.2 Arındırma sonrası (22 madde)

| Yöntem | Önerilen k |
|---|---|
| Velicer MAP | **1** |
| Paralel analiz — PC (%95) | **1** |
| Paralel analiz — indirgenmiş matris (%95) | 3 (marjinal: 1,173 vs. 1,048) |
| 1./2. özdeğer oranı | 7,17 (> 4 → tek boyutlu) |
| Bifaktör ω_H | 0,83 (> 0,70) |
| Bifaktör ECV | 0,749 (> 0,70) |
| ω_H / ω_total | 0,869 |

Rodriguez, Reise & Haviland (2016) ölçütlerine göre (ECV > 0,70 **ve** ω_H > 0,70) yapı
**özsel olarak tek boyutludur.**

---

## 5. Madde Çıkarma Süreci

Toplam **9 madde** iki ayrı gerekçeyle çıkarılmıştır.

### 5.1 Aşırı tavan etkisi (6 madde)
**Ölçüt:** Tepe kategoride ≥ %85 yanıt **ve** |çarpıklık| > 3 (Kline, 2016: |çarpıklık| > 3 ve
|basıklık| > 10 ciddi ihlal).

| Madde | Ort. | Tepe kat. % | Çarpıklık | Basıklık | Bootstrap SE |
|---|---|---|---|---|---|
| M14 | 4,82 | 91,6 | −4,45 | 19,75 | 0,113 |
| M18 | 4,82 | 91,3 | −4,44 | 20,32 | 0,121 |
| M20 | 4,77 | 89,9 | −3,57 | 12,30 | 0,112 |
| M4 | 4,74 | 88,7 | −3,55 | 11,93 | 0,102 |
| M3 | 4,72 | 87,5 | −3,18 | 9,41 | 0,099 |
| M12 | 4,82 | 86,8 | −4,06 | 19,66 | 0,103 |
| *(kalan maddeler)* | | | | | *0,081* |

Bu maddeler yalnızca dağılımsal olarak sorunlu değildir; somut zararları belgelenmiştir:

1. **Ayırt edicilik yok:** Katılımcıların ~%90'ı aynı kategoriyi işaretlemektedir.
2. **Kestirim istikrarsızlığı:** Bootstrap standart hataları diğer maddelerden %33 daha yüksektir;
   en gürültülü 10 korelasyon çiftinin tamamı bu maddeleri içermektedir. M18'in bir yanıt
   kategorisi tamamen **boştur** (n = 0).
3. **Tekilliğe yol açma:** 31 maddelik polikorik matris pozitif tanımlı değildir (3 negatif
   özdeğer); en negatif özdeğerin özvektörü M18 (0,58) ve M20 (0,37) tarafından domine edilmektedir.
4. **Sahte faktör üretme:** M18–M20 arasındaki şişmiş r = 0,887, k = 4 ve k = 5 çözümlerinde
   yalnızca bu iki maddeden oluşan yapay bir "ikili faktör" (yükler 0,86 ve 0,97) yaratmaktadır.
   Bu, çarpık maddelerin ürettiği klasik **güçlük faktörü** artefaktıdır (Bernstein & Teng, 1989).

Çıkarılmalarının ardından matris pozitif tanımlı hâle gelmiş, en yüksek korelasyon 0,887'den
0,650'ye düşmüş, KMO 0,852 → 0,882 yükselmiştir.

### 5.2 Yapıyla ilişkisizlik (3 madde)
**Ölçüt:** Diğer maddelerle ihmal edilebilir ilişki + düşük madde-toplam korelasyonu + düşük MSA.

| Madde | Diğer 20 maddeyle ort. r | Madde-toplam r | MSA | h² |
|---|---|---|---|---|
| M11 | **0,081** | 0,26 | 0,567 | 0,12–0,17 |
| M17 | **0,025** | 0,20 | 0,690 | 0,24–0,28 |
| M29 | **0,049** | 0,21 | 0,659 | 0,27–0,32 |
| *(diğer maddeler kendi aralarında)* | *0,449* | | | |

Bu üç madde diğerleriyle pratikte **hiç ilişkili değildir** (r ≈ 0,03–0,08; diğerleri kendi
aralarında 0,45). Kendi aralarındaki ilişkileri de zayıftır (0,20–0,30) ve ölçülebilir bir
faktör oluşturmaya yetmemektedir. Üçü de MSA < 0,70 ve madde-toplam r < 0,30 ölçütlerini
ihlal etmektedir (Kaiser & Rice, 1974; Field, 2013).

Çıkarılmalarının ardından KMO 0,882 → 0,899, en düşük MSA 0,246 → 0,664 olmuştur.

---

## 6. Final Çözüm: 22 Madde, Tek Faktör

**Çıkarım:** MINRES/ULS, polikorik korelasyon · **Açıklanan varyans:** %45,4
**Özdeğerler:** 10,51 · 1,47 · 1,26 · 0,97 · 0,89

| Madde | Ort. | SS | Faktör yükü | h² | MSA | Madde-toplam r | α (madde atılırsa) |
|---|---|---|---|---|---|---|---|
| M27 | 3,89 | 1,27 | **0,821** | 0,674 | 0,849 | 0,700 | 0,943 |
| M21 | 4,31 | 1,06 | **0,756** | 0,572 | 0,848 | 0,594 | 0,943 |
| M1 | 4,10 | 1,26 | **0,735** | 0,540 | 0,879 | 0,603 | 0,944 |
| M22 | 3,95 | 1,42 | **0,730** | 0,533 | 0,919 | 0,587 | 0,944 |
| M31 | 4,73 | 0,73 | **0,720** | 0,519 | 0,801 | 0,499 | 0,944 |
| M8 | 4,50 | 0,97 | **0,718** | 0,516 | 0,847 | 0,515 | 0,944 |
| M30 | 4,63 | 0,78 | **0,718** | 0,516 | 0,879 | 0,510 | 0,944 |
| M10 | 4,22 | 1,24 | **0,710** | 0,503 | 0,786 | 0,568 | 0,944 |
| M24 | 3,04 | 1,77 | **0,696** | 0,485 | 0,704 | 0,583 | 0,944 |
| M26 | 3,89 | 1,40 | **0,694** | 0,482 | 0,825 | 0,568 | 0,944 |
| M2 | 4,58 | 0,90 | **0,693** | 0,481 | 0,755 | 0,490 | 0,944 |
| M7 | 4,38 | 1,09 | **0,689** | 0,475 | 0,933 | 0,533 | 0,944 |
| M9ᴿ | 3,98 | 1,37 | **0,663** | 0,440 | 0,797 | 0,534 | 0,945 |
| M6 | 4,60 | 0,89 | **0,656** | 0,430 | 0,721 | 0,454 | 0,945 |
| M28 | 4,59 | 0,88 | **0,651** | 0,423 | 0,939 | 0,455 | 0,945 |
| M19 | 4,00 | 1,41 | **0,640** | 0,409 | 0,821 | 0,507 | 0,945 |
| M23 | 4,42 | 1,11 | **0,639** | 0,409 | 0,933 | 0,450 | 0,945 |
| M13 | 4,01 | 1,31 | **0,592** | 0,350 | 0,805 | 0,480 | 0,946 |
| M25 | 4,72 | 0,75 | **0,584** | 0,342 | 0,708 | 0,364 | 0,946 |
| M5 | 4,67 | 0,71 | **0,570** | 0,325 | 0,664 | 0,382 | 0,946 |
| M16 | 4,48 | 1,10 | **0,535** | 0,287 | 0,712 | 0,346 | 0,946 |
| M15 | 4,58 | 1,00 | **0,523** | 0,274 | 0,790 | 0,310 | 0,947 |

ᴿ = ters kodlanmış

- Tüm yükler ≥ 0,52 (Hair vd., 2019 ölçütü: ≥ 0,50 pratik anlamlılık)
- Hiçbir madde atıldığında α artmamaktadır → tüm maddeler ölçeğe katkı sağlamaktadır
- **Güvenirlik:** Sıralı α = 0,947 · McDonald ω = 0,948
- **Ölçek puanı:** Ort = 4,28 · SS = 0,64 · çarpıklık = −0,37 · basıklık = −0,81
  (toplam puan düzeyinde normallik sağlanmıştır)

### Uyum iyiliği (robust, WLSMV tipi)
```
χ²(mv) = 110,9   df = 64,8   p < 0,001
RMSEA = 0,048    CFI = 0,969    TLI = 0,982    SRMR = 0,077
```
Tüm indeksler kabul ölçütlerini karşılamaktadır (Hu & Bentler, 1999).

---

## 7. Neden Tek Faktör? — Aday Modellerin Karşılaştırması

| Model | Madde | k | RMSEA | CFI | TLI | SRMR | En düşük yük | Binişik | Φ_maks | **Kararlılık (φ)** |
|---|---|---|---|---|---|---|---|---|---|---|
| 31 madde, 4 faktör | 31 | 4 | 0,023 | 0,987 | 0,992 | 0,060 | 0,29 | 7 | 0,50 | — |
| 25 madde, 3 faktör | 25 | 3 | 0,034 | 0,982 | 0,988 | 0,059 | 0,41 | 9 | 0,48 | — |
| 22 madde, 3 faktör | 22 | 3 | 0,023 | 0,994 | 0,996 | 0,053 | **0,24** | 2 | **0,68** | **0,648** |
| 22 madde, 2 faktör | 22 | 2 | 0,038 | 0,982 | 0,989 | 0,064 | 0,37 | 1 | **0,70** | **0,753** |
| **22 madde, TEK faktör** | **22** | **1** | **0,048** | **0,969** | **0,982** | **0,077** | **0,52** | **0** | — | **0,993** |
| 20 madde, tek faktör | 20 | 1 | 0,053 | 0,967 | 0,981 | 0,076 | 0,57 | 0 | — | — |

**Yorum.** Uyum indeksleri faktör sayısıyla birlikte matematiksel zorunlulukla iyileşir; bu
nedenle tek başına uyum, boyutluluk kararı için yeterli değildir. Belirleyici kanıt
**çapraz geçerlemedir:**

| k | Ortalama Tucker φ | φ ≥ 0,85 oranı | En zayıf faktör φ |
|---|---|---|---|
| **1** | **0,993** | **%100** | 0,993 |
| 2 | 0,753 | %23 | 0,695 |
| 3 | 0,648 | %8 | 0,536 |
| 4 | 0,630 | %10 | 0,463 |

*(100 rastgele yarı-örneklem; Lorenzo-Seva & ten Berge, 2006: φ ≥ 0,95 eşdeğer, 0,85–0,94 benzer,
< 0,85 kararsız)*

Çok faktörlü çözümlerdeki ikinci ve üçüncü faktörler **örneklemden örnekleme yeniden
üretilememektedir.** Ayrıca k = 2 ve k = 3'te faktörler arası korelasyon 0,68–0,70 düzeyindedir
(ayırt edici geçerlik zayıf) ve k = 3'te M16 hiçbir faktöre 0,32 üzerinde yüklenmemektedir.
Buna karşılık tek faktörlü yapı neredeyse kusursuz replike olmaktadır.

---

## 8. Öneriler

1. **Ölçek tek boyutlu olarak puanlanmalıdır** (22 maddenin ortalaması veya toplamı).
   Alt boyut puanı hesaplanması önerilmez.
2. **M9'un ters kodlanmış olduğundan emin olunuz.**
3. α = 0,947 oldukça yüksektir; madde içeriklerinde gereksiz tekrar (redundans) olup olmadığı
   gözden geçirilebilir. Kısa form gerekiyorsa en yüksek yüklü 10 madde
   (M27, M21, M1, M22, M31, M8, M30, M10, M24, M26) iyi bir başlangıç noktasıdır.
4. Bulgular **açımlayıcıdır**; yapının **bağımsız bir örneklemde DFA** ile sınanması gerekir.
5. Atılan 6 tavan etkili maddenin içerik açısından kritik olup olmadığı değerlendirilmelidir.
   Kritikse, bu maddeler daha fazla değişkenlik üretecek biçimde (ör. daha geniş yanıt ölçeği,
   daha ayırt edici ifade) yeniden yazılabilir.

## 9. Sınırlılıklar

- **Madde içerikleri elimizde yoktur.** Tüm kararlar istatistikseldir; kapsam geçerliği ve
  içerik yorumu araştırmacı tarafından değerlendirilmelidir. Özellikle atılan maddelerin
  kuramsal olarak vazgeçilmez olup olmadığı içerik temelinde denetlenmelidir.
- n = 310, 22 madde için madde başına ~14 katılımcıya karşılık gelir; yeterlidir ancak
  çok faktörlü çözümlerin kararsızlığında örneklem büyüklüğünün payı olabilir.
- Robust uyum indeksleri, Γ matrisinin bootstrap kestirimine dayanan WLSMV *tipi* bir
  düzeltmedir; Mplus/lavaan'ın analitik WLSMV uygulamasıyla birebir aynı değildir (simülasyonla
  doğrulanmıştır, ancak küçük sapmalar olabilir).
- Aynı örneklemde yürütülen çapraz geçerleme, bağımsız örneklem replikasyonunun yerini tutmaz.

---

## 10. Yeniden Üretilebilirlik

| Dosya | İçerik |
|---|---|
| `polychor.py` | Polikorik korelasyon (iki aşamalı ML, Gauss-Legendre bivariate normal CDF) |
| `rotate.py` | GPA eğik/dik döndürme (Jennrich, 2002) |
| `efa_core.py` | Çıkarım, uyum indeksleri, MAP testi |
| `robust_fit.py` | WLSMV tipi robust uyum indeksleri |
| `final_engine.py` | Madde ayıklama ölçütleri |
| `05_parallel.py`, `21_pa_final.py` | Paralel analiz |
| `16_gamma.py` | Bootstrap Γ (1000 tekrar) |
| `24_stability.py` | Yarı-örneklem kararlılık analizi |
| `25_final_solution.py` | Final çözüm |
| `final_22_madde.csv` | Final madde istatistikleri |

Doğrulama betikleri: `07_validate.py` (uyum indeksleri), `12_validate_poly.py` (polikorik
tahminci), `19_validate_robust.py` (robust düzeltme).

---

## Kaynaklar

Asparouhov, T., & Muthén, B. (2010). *Simple second order chi-square correction*. Mplus Technical Appendix.
Bernstein, I. H., & Teng, G. (1989). Factoring items and factoring scales are different. *Psychological Bulletin, 105*(3), 467–477.
Braeken, J., & van Assen, M. A. L. M. (2017). An empirical Kaiser criterion. *Psychological Methods, 22*(3), 450–466.
Field, A. (2013). *Discovering statistics using IBM SPSS statistics* (4th ed.). Sage.
Flora, D. B., & Curran, P. J. (2004). An empirical evaluation of alternative methods of estimation for confirmatory factor analysis with ordinal data. *Psychological Methods, 9*(4), 466–491.
Hair, J. F., Black, W. C., Babin, B. J., & Anderson, R. E. (2019). *Multivariate data analysis* (8th ed.). Cengage.
Holgado-Tello, F. P., Chacón-Moscoso, S., Barbero-García, I., & Vila-Abad, E. (2010). Polychoric versus Pearson correlations in EFA and CFA of ordinal variables. *Quality & Quantity, 44*(1), 153–166.
Hu, L., & Bentler, P. M. (1999). Cutoff criteria for fit indexes in covariance structure analysis. *Structural Equation Modeling, 6*(1), 1–55.
Jennrich, R. I. (2002). A simple general method for oblique rotation. *Psychometrika, 67*(1), 7–19.
Kaiser, H. F., & Rice, J. (1974). Little Jiffy, Mark IV. *Educational and Psychological Measurement, 34*(1), 111–117.
Kline, R. B. (2016). *Principles and practice of structural equation modeling* (4th ed.). Guilford.
Lorenzo-Seva, U., & ten Berge, J. M. F. (2006). Tucker's congruence coefficient as a meaningful index of factor similarity. *Methodology, 2*(2), 57–64.
Muthén, B. (1984). A general structural equation model with dichotomous, ordered categorical, and continuous latent variable indicators. *Psychometrika, 49*(1), 115–132.
Olsson, U. (1979). Maximum likelihood estimation of the polychoric correlation coefficient. *Psychometrika, 44*(4), 443–460.
Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor models. *Psychological Methods, 21*(2), 137–150.
Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics and standard errors in covariance structure analysis. In A. von Eye & C. C. Clogg (Eds.), *Latent variables analysis* (pp. 399–419). Sage.
Velicer, W. F. (1976). Determining the number of components from the matrix of partial correlations. *Psychometrika, 41*(3), 321–327.
