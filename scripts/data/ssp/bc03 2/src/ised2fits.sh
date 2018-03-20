#!/bin/sh

set name = $1:s|.ised||
$bc03/galaxevpl $name -all << !


!
$STILTS/stilts.sh   -verbose  tcopy  ifmt=ascii   ofmt=fits   $name.all   $name.fits
\rm $name.all
