# Veri sutunu -> yeni madde no (kontrol maddesi cikarilmis, 24-32 -> 23-31)
ITEM_NO = {('M%d'%i): i for i in range(1,32)}
ORIG_NO = {('M%d'%i): (i if i<=22 else i+1) for i in range(1,32)}
TXT = {
1:"Okul sorunlarini konusurken rahat hissederim",
2:"Sorun cozme yollarini bulmama yardim eder",
3:"[T] Basarili olup olmadigimi umursamaz",
4:"[T] Ogrenme surecimi umursamaz",
5:"Basarili olmam icin yuksek beklentileri var",
6:"Universiteye gidecegime dair yuksek beklenti",
7:"Sinavda basarili olunca over",
8:"Basarisiz olunca daha cok calismaya tesvik eder",
9:"[T] Basarili olmam icin buyuk baski hissederim",
10:"Gelecekteki okul planlarimi konusur",
11:"Meslek planlarimi konusur",
12:"Elimden gelenin en iyisini yapmami ister",
13:"Davranislarimi etkiler",
14:"[T] Okula gidip gitmedigimle ilgilenmez",
15:"[T] Ogretmenlerimle hic gorusmedi",
16:"Yardimci kaynak icin para verir",
17:"Ek ogrenme materyalleri alir",
18:"[T] Ek materyal karsilayacak parasi yok",
19:"Ozel ders almam icin tesvik eder",
20:"[T] Ozel ders aldiracak parasi yok",
21:"Sik sik okul calismalarimi sorar",
22:"[T] Odevlerimde yardim edecek zamani yok",
23:"[T] Yardim edecek yeterli bilgisi yok",
24:"[T] Dusuk basarili arkadaslarla sosyallesmemi tavsiye eder",
25:"Iyi huylu arkadaslarla sosyallesmemi tavsiye eder",
26:"Ders zamanim olsun diye ev isi istemez",
27:"Zamanimin cogunu ders calisarak gecirdigimden emin olur",
28:"Karnelerimi dikkatle inceler",
29:"Sik sik notlarimi sorar",
30:"Odevimi nasil yapacagimi bildigime guvenir",
31:"GENEL: Basarili olmama yardim icin destek olur",
}
REVERSE = [3,4,9,14,15,18,20,22,23,24]      # yeni numaralandirmada ters maddeler
def lab(c, w=52):
    n=ITEM_NO[c]; return "%-4s m%-2d %s"%(c,n,TXT[n][:w])
