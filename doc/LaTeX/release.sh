#!/bin/sh
if [ -f document.pdf ]; then
  cp document.pdf document-$(date "+%Y-%m-%d").pdf
fi
