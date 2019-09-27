filename="$1"
while read line; do Rscript /well/mccarthy/users/maxlouis/oxford2/HNF1A_project/scripts/table2/2.count_1read.R $line; done < "$filename"
