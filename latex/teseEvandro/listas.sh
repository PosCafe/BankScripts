#!/bin/sh
makeindex -s teseEvandro.ist -t teseEvandro.sgl -o teseEvandro.sbl teseEvandro.sim
makeindex -s teseEvandro.ist -t teseEvandro.agl -o teseEvandro.sgl teseEvandro.sig


#find ./ -maxdepth 1 -name '*.ps' -exec pstoimg -type png -crop a -density 600 -aaliastext -color 24 -debug  {} \;
