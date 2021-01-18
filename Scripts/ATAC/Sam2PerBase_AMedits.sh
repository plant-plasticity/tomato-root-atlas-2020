#!/bin/bash 

awk -F"\t" '{if($3 == "Sl3.0ch01" || $3 == "SL3.0ch01") print $0 }' $1 | gawk '(! and(16, $2))' | awk '{print "SL3.0ch01\t"$4-1"\t"$4"\tSL3.0ch01"}' > chr1.Fperbase.bed
awk -F"\t" '{if($3 == "Sl3.0ch01" || $3 == "SL3.0ch01") print $0 }' $1 | gawk '(and(16, $2))' | awk '{print "SL3.0ch01\t"$4+35"\t"$4+36"\tSL3.0ch01"}' > chr1.Rperbase.bed
cat chr1.Fperbase.bed chr1.Rperbase.bed > chr1.FRperbase.bed
sort-bed chr1.FRperbase.bed > chr1.sorted_FRperbase.bed
bedmap --echo --bases-uniq --delim "\t" chr1.sorted_FRperbase.bed chr1.sorted_FRperbase.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' | uniq > chr1.merged.bed
bedmap --echo --bases --delim "\t" chr1.merged.bed chr1.sorted_FRperbase.bed > chr1.perbase.bed
rm chr1.Fperbase.bed
rm chr1.Rperbase.bed
rm chr1.FRperbase.bed
rm chr1.sorted_FRperbase.bed
rm chr1.merged.bed

awk -F"\t" '{if($3 == "Sl3.0ch02" || $3 == "SL3.0ch02") print $0 }' $1 | gawk '(! and(16, $2))' | awk '{print "SL3.0ch02\t"$4-1"\t"$4"\tSL3.0ch02"}' > chr2.Fperbase.bed
awk -F"\t" '{if($3 == "Sl3.0ch02" || $3 == "SL3.0ch02") print $0 }' $1 | gawk '(and(16, $2))' | awk '{print "SL3.0ch02\t"$4+35"\t"$4+36"\tSL3.0ch02"}' > chr2.Rperbase.bed
cat chr2.Fperbase.bed chr2.Rperbase.bed > chr2.FRperbase.bed
sort-bed chr2.FRperbase.bed > chr2.sorted_FRperbase.bed
bedmap --echo --bases-uniq --delim "\t" chr2.sorted_FRperbase.bed chr2.sorted_FRperbase.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' | uniq > chr2.merged.bed
bedmap --echo --bases --delim "\t" chr2.merged.bed chr2.sorted_FRperbase.bed > chr2.perbase.bed
rm chr2.Fperbase.bed
rm chr2.Rperbase.bed
rm chr2.FRperbase.bed
rm chr2.sorted_FRperbase.bed
rm chr2.merged.bed

awk -F"\t" '{if($3 == "Sl3.0ch03" || $3 == "SL3.0ch03") print $0 }' $1 | gawk '(! and(16, $2))' | awk '{print "SL3.0ch03\t"$4-1"\t"$4"\tSL3.0ch03"}' > chr3.Fperbase.bed
awk -F"\t" '{if($3 == "Sl3.0ch03" || $3 == "SL3.0ch03") print $0 }' $1 | gawk '(and(16, $2))' | awk '{print "SL3.0ch03\t"$4+35"\t"$4+36"\tSL3.0ch03"}' > chr3.Rperbase.bed
cat chr3.Fperbase.bed chr3.Rperbase.bed > chr3.FRperbase.bed
sort-bed chr3.FRperbase.bed > chr3.sorted_FRperbase.bed
bedmap --echo --bases-uniq --delim "\t" chr3.sorted_FRperbase.bed chr3.sorted_FRperbase.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' | uniq > chr3.merged.bed
bedmap --echo --bases --delim "\t" chr3.merged.bed chr3.sorted_FRperbase.bed > chr3.perbase.bed
rm chr3.Fperbase.bed
rm chr3.Rperbase.bed
rm chr3.FRperbase.bed
rm chr3.sorted_FRperbase.bed
rm chr3.merged.bed

awk -F"\t" '{if($3 == "Sl3.0ch04" || $3 == "SL3.0ch04") print $0 }' $1 | gawk '(! and(16, $2))' | awk '{print "SL3.0ch04\t"$4-1"\t"$4"\tSL3.0ch04"}' > chr4.Fperbase.bed
awk -F"\t" '{if($3 == "Sl3.0ch04" || $3 == "SL3.0ch04") print $0 }' $1 | gawk '(and(16, $2))' | awk '{print "SL3.0ch04\t"$4+35"\t"$4+36"\tSL3.0ch04"}' > chr4.Rperbase.bed
cat chr4.Fperbase.bed chr4.Rperbase.bed > chr4.FRperbase.bed
sort-bed chr4.FRperbase.bed > chr4.sorted_FRperbase.bed
bedmap --echo --bases-uniq --delim "\t" chr4.sorted_FRperbase.bed chr4.sorted_FRperbase.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' | uniq > chr4.merged.bed
bedmap --echo --bases --delim "\t" chr4.merged.bed chr4.sorted_FRperbase.bed > chr4.perbase.bed
rm chr4.Fperbase.bed
rm chr4.Rperbase.bed
rm chr4.FRperbase.bed
rm chr4.sorted_FRperbase.bed
rm chr4.merged.bed

awk -F"\t" '{if($3 == "Sl3.0ch05" || $3 == "SL3.0ch05") print $0 }' $1 | gawk '(! and(16, $2))' | awk '{print "SL3.0ch05\t"$4-1"\t"$4"\tSL3.0ch05"}' > chr5.Fperbase.bed
awk -F"\t" '{if($3 == "Sl3.0ch05" || $3 == "SL3.0ch05") print $0 }' $1 | gawk '(and(16, $2))' | awk '{print "SL3.0ch05\t"$4+35"\t"$4+36"\tSL3.0ch05"}' > chr5.Rperbase.bed
cat chr5.Fperbase.bed chr5.Rperbase.bed > chr5.FRperbase.bed
sort-bed chr5.FRperbase.bed > chr5.sorted_FRperbase.bed
bedmap --echo --bases-uniq --delim "\t" chr5.sorted_FRperbase.bed chr5.sorted_FRperbase.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' | uniq > chr5.merged.bed
bedmap --echo --bases --delim "\t" chr5.merged.bed chr5.sorted_FRperbase.bed > chr5.perbase.bed
rm chr5.Fperbase.bed
rm chr5.Rperbase.bed
rm chr5.FRperbase.bed
rm chr5.sorted_FRperbase.bed
rm chr5.merged.bed

awk -F"\t" '{if($3 == "Sl3.0ch06" || $3 == "SL3.0ch06") print $0 }' $1 | gawk '(! and(16, $2))' | awk '{print "SL3.0ch06\t"$4-1"\t"$4"\tSL3.0ch06"}' > chr6.Fperbase.bed
awk -F"\t" '{if($3 == "Sl3.0ch06" || $3 == "SL3.0ch06") print $0 }' $1 | gawk '(and(16, $2))' | awk '{print "SL3.0ch06\t"$4+35"\t"$4+36"\tSL3.0ch06"}' > chr6.Rperbase.bed
cat chr6.Fperbase.bed chr6.Rperbase.bed > chr6.FRperbase.bed
sort-bed chr6.FRperbase.bed > chr6.sorted_FRperbase.bed
bedmap --echo --bases-uniq --delim "\t" chr6.sorted_FRperbase.bed chr6.sorted_FRperbase.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' | uniq > chr6.merged.bed
bedmap --echo --bases --delim "\t" chr6.merged.bed chr6.sorted_FRperbase.bed > chr6.perbase.bed
rm chr6.Fperbase.bed
rm chr6.Rperbase.bed
rm chr6.FRperbase.bed
rm chr6.sorted_FRperbase.bed
rm chr6.merged.bed

awk -F"\t" '{if($3 == "Sl3.0ch07" || $3 == "SL3.0ch07") print $0 }' $1 | gawk '(! and(16, $2))' | awk '{print "SL3.0ch07\t"$4-1"\t"$4"\tSL3.0ch07"}' > chr7.Fperbase.bed
awk -F"\t" '{if($3 == "Sl3.0ch07" || $3 == "SL3.0ch07") print $0 }' $1 | gawk '(and(16, $2))' | awk '{print "SL3.0ch07\t"$4+35"\t"$4+36"\tSL3.0ch07"}' > chr7.Rperbase.bed
cat chr7.Fperbase.bed chr7.Rperbase.bed > chr7.FRperbase.bed
sort-bed chr7.FRperbase.bed > chr7.sorted_FRperbase.bed
bedmap --echo --bases-uniq --delim "\t" chr7.sorted_FRperbase.bed chr7.sorted_FRperbase.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' | uniq > chr7.merged.bed
bedmap --echo --bases --delim "\t" chr7.merged.bed chr7.sorted_FRperbase.bed > chr7.perbase.bed
rm chr7.Fperbase.bed
rm chr7.Rperbase.bed
rm chr7.FRperbase.bed
rm chr7.sorted_FRperbase.bed
rm chr7.merged.bed

awk -F"\t" '{if($3 == "Sl3.0ch08" || $3 == "SL3.0ch08") print $0 }' $1 | gawk '(! and(16, $2))' | awk '{print "SL3.0ch08\t"$4-1"\t"$4"\tSL3.0ch08"}' > chr8.Fperbase.bed
awk -F"\t" '{if($3 == "Sl3.0ch08" || $3 == "SL3.0ch08") print $0 }' $1 | gawk '(and(16, $2))' | awk '{print "SL3.0ch08\t"$4+35"\t"$4+36"\tSL3.0ch08"}' > chr8.Rperbase.bed
cat chr8.Fperbase.bed chr8.Rperbase.bed > chr8.FRperbase.bed
sort-bed chr8.FRperbase.bed > chr8.sorted_FRperbase.bed
bedmap --echo --bases-uniq --delim "\t" chr8.sorted_FRperbase.bed chr8.sorted_FRperbase.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' | uniq > chr8.merged.bed
bedmap --echo --bases --delim "\t" chr8.merged.bed chr8.sorted_FRperbase.bed > chr8.perbase.bed
rm chr8.Fperbase.bed
rm chr8.Rperbase.bed
rm chr8.FRperbase.bed
rm chr8.sorted_FRperbase.bed
rm chr8.merged.bed

awk -F"\t" '{if($3 == "Sl3.0ch09" || $3 == "SL3.0ch09") print $0 }' $1 | gawk '(! and(16, $2))' | awk '{print "SL3.0ch09\t"$4-1"\t"$4"\tSL3.0ch09"}' > chr9.Fperbase.bed
awk -F"\t" '{if($3 == "Sl3.0ch09" || $3 == "SL3.0ch09") print $0 }' $1 | gawk '(and(16, $2))' | awk '{print "SL3.0ch09\t"$4+35"\t"$4+36"\tSL3.0ch09"}' > chr9.Rperbase.bed
cat chr9.Fperbase.bed chr9.Rperbase.bed > chr9.FRperbase.bed
sort-bed chr9.FRperbase.bed > chr9.sorted_FRperbase.bed
bedmap --echo --bases-uniq --delim "\t" chr9.sorted_FRperbase.bed chr9.sorted_FRperbase.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' | uniq > chr9.merged.bed
bedmap --echo --bases --delim "\t" chr9.merged.bed chr9.sorted_FRperbase.bed > chr9.perbase.bed
rm chr9.Fperbase.bed
rm chr9.Rperbase.bed
rm chr9.FRperbase.bed
rm chr9.sorted_FRperbase.bed
rm chr9.merged.bed

awk -F"\t" '{if($3 == "Sl3.0ch10" || $3 == "SL3.0ch10") print $0 }' $1 | gawk '(! and(16, $2))' | awk '{print "SL3.0ch10\t"$4-1"\t"$4"\tSL3.0ch10"}' > chr10.Fperbase.bed
awk -F"\t" '{if($3 == "Sl3.0ch10" || $3 == "SL3.0ch10") print $0 }' $1 | gawk '(and(16, $2))' | awk '{print "SL3.0ch10\t"$4+35"\t"$4+36"\tSL3.0ch10"}' > chr10.Rperbase.bed
cat chr10.Fperbase.bed chr10.Rperbase.bed > chr10.FRperbase.bed
sort-bed chr10.FRperbase.bed > chr10.sorted_FRperbase.bed
bedmap --echo --bases-uniq --delim "\t" chr10.sorted_FRperbase.bed chr10.sorted_FRperbase.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' | uniq > chr10.merged.bed
bedmap --echo --bases --delim "\t" chr10.merged.bed chr10.sorted_FRperbase.bed > chr10.perbase.bed
rm chr10.Fperbase.bed
rm chr10.Rperbase.bed
rm chr10.FRperbase.bed
rm chr10.sorted_FRperbase.bed
rm chr10.merged.bed

awk -F"\t" '{if($3 == "Sl3.0ch11" || $3 == "SL3.0ch11") print $0 }' $1 | gawk '(! and(16, $2))' | awk '{print "SL3.0ch11\t"$4-1"\t"$4"\tSL3.0ch11"}' > chr11.Fperbase.bed
awk -F"\t" '{if($3 == "Sl3.0ch11" || $3 == "SL3.0ch11") print $0 }' $1 | gawk '(and(16, $2))' | awk '{print "SL3.0ch11\t"$4+35"\t"$4+36"\tSL3.0ch11"}' > chr11.Rperbase.bed
cat chr11.Fperbase.bed chr11.Rperbase.bed > chr11.FRperbase.bed
sort-bed chr11.FRperbase.bed > chr11.sorted_FRperbase.bed
bedmap --echo --bases-uniq --delim "\t" chr11.sorted_FRperbase.bed chr11.sorted_FRperbase.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' | uniq > chr11.merged.bed
bedmap --echo --bases --delim "\t" chr11.merged.bed chr11.sorted_FRperbase.bed > chr11.perbase.bed
rm chr11.Fperbase.bed
rm chr11.Rperbase.bed
rm chr11.FRperbase.bed
rm chr11.sorted_FRperbase.bed
rm chr11.merged.bed

awk -F"\t" '{if($3 == "Sl3.0ch12" || $3 == "SL3.0ch12") print $0 }' $1 | gawk '(! and(16, $2))' | awk '{print "SL3.0ch12\t"$4-1"\t"$4"\tSL3.0ch12"}' > chr12.Fperbase.bed
awk -F"\t" '{if($3 == "Sl3.0ch12" || $3 == "SL3.0ch12") print $0 }' $1 | gawk '(and(16, $2))' | awk '{print "SL3.0ch12\t"$4+35"\t"$4+36"\tSL3.0ch12"}' > chr12.Rperbase.bed
cat chr12.Fperbase.bed chr12.Rperbase.bed > chr12.FRperbase.bed
sort-bed chr12.FRperbase.bed > chr12.sorted_FRperbase.bed
bedmap --echo --bases-uniq --delim "\t" chr12.sorted_FRperbase.bed chr12.sorted_FRperbase.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' | uniq > chr12.merged.bed
bedmap --echo --bases --delim "\t" chr12.merged.bed chr12.sorted_FRperbase.bed > chr12.perbase.bed
rm chr12.Fperbase.bed
rm chr12.Rperbase.bed
rm chr12.FRperbase.bed
rm chr12.sorted_FRperbase.bed
rm chr12.merged.bed

sort-bed chr1.perbase.bed > chr1.sorted_perbase.bed
sort-bed chr2.perbase.bed > chr2.sorted_perbase.bed
sort-bed chr3.perbase.bed > chr3.sorted_perbase.bed
sort-bed chr4.perbase.bed > chr4.sorted_perbase.bed
sort-bed chr5.perbase.bed > chr5.sorted_perbase.bed
sort-bed chr6.perbase.bed > chr6.sorted_perbase.bed
sort-bed chr7.perbase.bed > chr7.sorted_perbase.bed
sort-bed chr8.perbase.bed > chr8.sorted_perbase.bed
sort-bed chr9.perbase.bed > chr9.sorted_perbase.bed
sort-bed chr10.perbase.bed > chr10.sorted_perbase.bed
sort-bed chr11.perbase.bed > chr11.sorted_perbase.bed
sort-bed chr12.perbase.bed > chr12.sorted_perbase.bed

cat chr1.sorted_perbase.bed chr2.sorted_perbase.bed chr3.sorted_perbase.bed chr4.sorted_perbase.bed chr5.sorted_perbase.bed chr6.sorted_perbase.bed chr7.sorted_perbase.bed chr8.sorted_perbase.bed chr9.sorted_perbase.bed chr10.sorted_perbase.bed chr11.sorted_perbase.bed chr12.sorted_perbase.bed > big.perbase

rm chr1.perbase.bed
rm chr2.perbase.bed
rm chr3.perbase.bed
rm chr4.perbase.bed
rm chr5.perbase.bed
rm chr6.perbase.bed
rm chr7.perbase.bed
rm chr8.perbase.bed
rm chr9.perbase.bed
rm chr10.perbase.bed
rm chr11.perbase.bed
rm chr12.perbase.bed

rm chr1.sorted_perbase.bed
rm chr2.sorted_perbase.bed
rm chr3.sorted_perbase.bed
rm chr4.sorted_perbase.bed
rm chr5.sorted_perbase.bed
rm chr6.sorted_perbase.bed
rm chr7.sorted_perbase.bed
rm chr8.sorted_perbase.bed
rm chr9.sorted_perbase.bed
rm chr10.sorted_perbase.bed
rm chr11.sorted_perbase.bed
rm chr12.sorted_perbase.bed
