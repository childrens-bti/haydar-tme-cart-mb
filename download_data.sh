set -e
set -o pipefail

# Define URL and version
URL=${URL:-https://bti-openaccess-us-east-1-bti-bfx.s3.us-east-1.amazonaws.com/haydar-scrna}
RELEASE=${RELEASE:-v4}

# Remove old symlinks in data
find data -type l -delete

# Ensure release dir exists
mkdir -p "data/$RELEASE"

# The md5sum file provides our single point of truth for which files are in a release.
MD5="data/$RELEASE/md5sum.txt"
if [ -f "$MD5" ]; then
  curl -k --create-dirs -z "$MD5" -o "$MD5" "$URL/$RELEASE/md5sum.txt"
else
  # first run: no -z (no conditional GET)
  curl -k --create-dirs -o "$MD5" "$URL/$RELEASE/md5sum.txt"
fi

# Consider the filenames in the md5sum file
FILES=($(tr -s ' ' < "data/$RELEASE/md5sum.txt" \
         | cut -d ' ' -f 2 \
         | sed -e 's/^\*//' -e 's#^\./##'))

# Before each download, ensure parent directories exist
for file in "${FILES[@]}"; do
  if [ ! -e "data/$RELEASE/$file" ]; then
    mkdir -p "data/$RELEASE/$(dirname "$file")"
    echo "Downloading $file"
    curl -k "$URL/$RELEASE/$file" -o "data/$RELEASE/$file"
  fi
done

# Check the md5s for everything we downloaded except CHANGELOG.md
cd data/$RELEASE
echo "Checking MD5 hashes..."
md5sum -c md5sum.txt
cd ../../

# Make symlinks in data/ to the files in the just downloaded release folder.
#for file in "${FILES[@]}"
#do
#  ln -sfn $RELEASE/$file data/$file
#done