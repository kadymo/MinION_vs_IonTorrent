#1st line = gives you nr of sequences in the file
#2nd line = change # in {} , this then renames each sequence as e.g. AB1_seq# 


sed '/^>/d' AB1_19_L001_R1_001.fna | wc -l
for i in {1..111835}; do echo AB1_$i; done | paste - <(sed '/^>/d' AB1_19_L001_R1_001.fna) | sed -e 's/^/>/' -e 's/\t/\n/' > AB1.fna

sed '/^>/d' DNAs1.fna | wc -l
for i in {1..7872}; do echo AB1_$i; done | paste - <(sed '/^>/d' DNAs1.fna) | sed -e 's/^/>/' -e 's/\t/\n/' > DNAs1_new.fna

sed '/^>/d' DNAs2.fna | wc -l
for i in {1..7872}; do echo AB1_$i; done | paste - <(sed '/^>/d' DNAs2.fna) | sed -e 's/^/>/' -e 's/\t/\n/' > DNAs2_new.fna

sed '/^>/d' N42.fna | wc -l
for i in {1..116405}; do echo AB1_$i; done | paste - <(sed '/^>/d' N42.fna) | sed -e 's/^/>/' -e 's/\t/\n/' > N42_new.fna

sed '/^>/d' ODK3b.fna | wc -l
for i in {1..37772}; do echo AB1_$i; done | paste - <(sed '/^>/d' ODK3b.fna) | sed -e 's/^/>/' -e 's/\t/\n/' > ODK3b_new.fna

sed '/^>/d' ODK5b.fna | wc -l
for i in {1..25275}; do echo AB1_$i; done | paste - <(sed '/^>/d' ODK5b.fna) | sed -e 's/^/>/' -e 's/\t/\n/' > ODK5b_new.fna

sed '/^>/d' ODK7b.fna | wc -l
for i in {1..29106}; do echo AB1_$i; done | paste - <(sed '/^>/d' ODK7b.fna) | sed -e 's/^/>/' -e 's/\t/\n/' > ODK7b_new.fna

sed '/^>/d' ODK9a.fna | wc -l
for i in {1..36520}; do echo AB1_$i; done | paste - <(sed '/^>/d' ODK9a.fna) | sed -e 's/^/>/' -e 's/\t/\n/' > ODK9a_new.fna

sed '/^>/d' ODK11b.fna | wc -l
for i in {1..25073}; do echo AB1_$i; done | paste - <(sed '/^>/d' ODK11b.fna) | sed -e 's/^/>/' -e 's/\t/\n/' > ODK11b_new.fna

sed '/^>/d' ODK13b.fna | wc -l
for i in {1..40626}; do echo AB1_$i; done | paste - <(sed '/^>/d' ODK13b.fna) | sed -e 's/^/>/' -e 's/\t/\n/' > ODK13b_new.fna

sed '/^>/d' ODK15a.fna | wc -l
for i in {1..17923}; do echo AB1_$i; done | paste - <(sed '/^>/d' ODK15a.fna) | sed -e 's/^/>/' -e 's/\t/\n/' > ODK15a_new.fna

sed '/^>/d' ODK21b.fna | wc -l
for i in {1..31609}; do echo AB1_$i; done | paste - <(sed '/^>/d' ODK21b.fna) | sed -e 's/^/>/' -e 's/\t/\n/' > ODK21b_new.fna

sed '/^>/d' ODK36a.fna | wc -l
for i in {1..21563}; do echo AB1_$i; done | paste - <(sed '/^>/d' ODK36a.fna) | sed -e 's/^/>/' -e 's/\t/\n/' > ODK36a_new.fna
