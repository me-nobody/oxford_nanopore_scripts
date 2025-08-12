# the .detect format is huge ~GBs and hence i am using bash shell scripts to analyze it

# EBT
# count the total number of postions
grep -v "#" ebt_u2os_15112023_CPU.detect|wc -l
# 117881605
# EdU positive sites
grep -v "#" ebt_u2os_15112023_CPU.detect|awk -F "\t" '$2>0.5 {print $0}'|wc -l
# 175640
# % (175640/117881605) * 100 = 0.14 %
# BrdU positive sites
grep -v "#" ebt_u2os_15112023_CPU.detect|awk -F "\t" '$3>0.5 {print $0}'|wc -l
# 175869
# % (175869/117881605) * 100 = 0.14 %

# positive control
# count the total number of postions
grep -v "#" posctl_271123_CPU.detect|wc -l
# 103364352
# EdU positive sites
grep -v "#" posctl_271123_CPU.detect|awk -F "\t" '$2>0.5 {print $0}'|wc -l
# 655709

# BrdU positive sites
grep -v "#" posctl_271123_CPU.detect|awk -F "\t" '$3>0.5 {print $0}'|wc -l
# 36971
