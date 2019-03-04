#!/bin/sh
NS="$1"
shift
NS=`echo "$NS" | tr '[A-Z]' '[a-z]'`
NS1=`echo "$NS" | sed 's/\b./\u&/g'`
for p in "$@"; do
  cat > "${NS1}:${p}.txt" <<EOF
[[Imported from::${NS}:${p}]] 
{{DISPLAYTITLE:${NS}:${p}}}
[[Category:Imported vocabulary]]
EOF
done
