
{
  octave --quiet tests/UtilTest/conv/rand3.m | tee tmpP1.dat
  octave --quiet tests/UtilTest/conv/rand3.m | tee tmpP2.dat
} &> /dev/null

octave --quiet tests/UtilTest/conv/convOct.m
