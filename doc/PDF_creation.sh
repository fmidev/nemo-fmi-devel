#!/bin/sh


export opts='-shell-escape -interaction=nonstopmode'
model='NEMO'

clean() {
    ## Delete latex build files
    find latex -regextype posix-extended -regex ".*\.(aux|log|maf|mtc|out|toc).*" -exec rm {} \;

    ## Remove 'minted' directories
    find latex -type d -name '_minted*' -exec rm -r {} \;

    ## HTML exports
    find latex -type d -name 'html*'    -exec rm -r {} \;
}

build() {
    cd latex/$1/main

    latex     $opts              $1'_manual' > /dev/null
    makeindex -s $1'_manual'.ist $1'_manual' > /dev/null
    bibtex                       $1'_manual' > /dev/null
    #latex     $opts              $1'_manual' > /dev/null
    pdflatex  $opts              $1'_manual' > /dev/null

    mv $1'_manual'.pdf ../../..
    cd -
}

clean
build    $model

exit 0
