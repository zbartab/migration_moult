R --vanilla --no-save <<EOF
library(Rmoult,lib.loc="~/lib/R/library");
postscript(file="$1.ps",horizontal=F);
sink(file="$1.txt")
p.ser(patt="^$1.*\\.dat\$",ask=F);
sink();
dev.off();
EOF