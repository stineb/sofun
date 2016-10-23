#!/bin/sh
#
# batch_audio_cvt.sh
#
# written by Tyler W. Davis
#
# 2015-07-04 -- created
# 2016-10-23 -- last updated
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script performs batch audio conversion from m4a to mp3.
#
# ~~~~~~~~~~~~~
# requirements:
# ~~~~~~~~~~~~~
# libav-tools (e.g., v.6:9.18-0ubuntu0.14.04)
# mediainfo (e.g, v.0.7.67)
#
# ~~~~~~~~~~~
# references:
# ~~~~~~~~~~~
# http://giantdorks.org/alain/a-few-methods-to-convert-aac-m4a-to-mp3-under-linux/
#
# ~~~~~~
# usage:
# ~~~~~~
# `bash batch_audio_cvt.sh *.m4a`

GetInfo()
{
  info=$(mediainfo "$f" | sed -r 's/[ ]+:[ ]/|/') 
  bitrate=$(echo "$info" | awk -F'|' '$1=="Bit rate" {gsub(/([ A-Za-z]+)/,"",$2); print$2}')
}

for f in "$@"; do
  GetInfo
  avconv -i "$f" -ab ${bitrate}k "${f%.m4a}.mp3"
done

