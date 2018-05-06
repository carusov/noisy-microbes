#!/usr/bin/awk -f

BEGIN {
    OFS = "\t"
}

/^>/ {
    sub(">", "", $1)
#    gsub(":", "_", $1)
#    sub("_[0-9]*", "", $1)
    n=split(FILENAME, path, "/")
    file=path[n]
    split(file, sample, "_")
#    gsub("-", "_", sample[1])
    print $1, sample[1]
}
