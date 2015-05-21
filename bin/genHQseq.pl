$i = 1;

open MATURE, '>', "hq_matureMirna.fas";
open HAIRPIN, '>', "hq_hairpinMirna.fas";

while(<>){
    chomp;
    ($uniqName, $knownName, $matureSeq, $hairpinSeq) = split "\t", $_;
    if( $knownName eq "-" ){
        print MATURE ">cro-mir-T$i\n$matureSeq\n";
        print HAIRPIN ">cro-mir-T$i\n$hairpinSeq\n";
    }else{
        print MATURE ">cro-mir-K$i\n$matureSeq\n";
        print HAIRPIN ">cro-mir-K$i\n$hairpinSeq\n";
    }
    $i++;
}
