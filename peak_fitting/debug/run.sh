set -e
trap 'case $? in
        139) echo "segfault occurred";;
      esac' EXIT
for i in {1..100}
  do
  Rscript SA.wrapper.r 19
  done

