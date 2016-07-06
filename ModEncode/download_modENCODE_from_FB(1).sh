#!/bin/sh

#./download_modENCODE_from_FB.sh dmel_FBgn_sorted dmel_modENCODE_raw.tsv

while read line
do
  wget wget http\:\/\/flybase.org\/cgi-bin\/serveHTdata\.cgi\?dataset\=modENCODE_mRNA-Seq_tissues\&FBgn\="$line" -O - >> "$2".2
done < $1

#headers
grep -v "^#" "$2".2 > "$2"
head -n1 "$2".2 | perl -pe 's/#//g' | perl -pe 's/ /_/g' | cat - "$2" > "$2".3
mv "$2".3 "$2"
rm "$2".2

