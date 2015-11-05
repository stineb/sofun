#!/bin/sh
#
# batch_audio_cvt.sh
#
# written by Tyler W. Davis
#
# 2015-07-04 -- created
# 2015-07-04 -- last updated
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script performs batch audio/video conversion based on VLC Media Player
# (https://www.videolan.org/vlc/). For example, to batch convert lossless m4a 
# audio to lossly mp3: `bash batch_audio_cvt.sh *.m4a`
#
# ~~~~~~~~~~~
# references:
# ~~~~~~~~~~~
# "How to Batch Encode," online: https://wiki.videolan.org/How_to_Batch_Encode
#
mux="mp3"
vlc="/usr/bin/vlc"

if [ ! -e "$vlc" ]; then
    echo "Command '$vlc' does not exist"
    exit 1
fi

for file in "$@"; do
    echo "=> Transcoding '$file'... "

    # NOTE: basename is setup for letter-number-letter format (e.g., m4a)
    dst=`dirname "$file"`
    new=`basename "$file" | sed 's@\.[a-z][0-9][a-z]$@@'`.$mux

    $vlc -I dummy -q "$file" \
       --sout "#transcode{acodec=mp4a,ab=128}:standard{mux=mp4,dst=\"$dst/$new\",access=file}" \
       vlc://quit
    ls -lh "$file" "$dst/$new"
    echo
done
