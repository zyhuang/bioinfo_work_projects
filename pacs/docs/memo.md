### Notes on Annovar Gene-based Annotations

[Reference](http://annovar.openbioinformatics.org/en/latest/user-guide/gene/)

`Func.refGene` values in order of descending precedence (similar for `Func.wgEncodeGencodeBasicV19`). In case of "both", the values are separated by `\x3b` (`;`)

* `exonic` / `splicing` / both (=> see `ExonicFunc.refGene`)
* `ncRNA_exonic` / `ncRNA_splicing` / both 
* `UTR5` / `UTR3` / both 
* `intronic`
* `ncRNA_intronic`
* `upstream` / `downstream` / both 
* `intergenic`

When `Func.refGene` contains `exonic`, `ExonicFunc.refGene` values in order of descending precedence (similar for `ExonicFunc.wgEncodeGencodeBasicV19`)

* `frameshift_insertion`
* `frameshift_deletion`
* `frameshift_block_substitution`
* `stopgain`
* `stoploss`
* `nonframeshift_insertion`
* `nonframeshift_deletion`
* `nonframeshift_block_substitution`
* `nonsynonymous_SNV`
* `synonymous_SNV`
* `unknown`



