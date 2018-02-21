This file documents the format of intermediate data in PACS.

### variant calling output 

* data location: `data2/pacs/var_json/batch_[BID]/[SID].var.json.gz` (Batch ID `BID` = `00..21`, `SID` is sample ID).
* data size: 4.7TB. 
* number of files: 214277 (1 per sample, 10000 per batch).
* content: Random Forest variant classification result, one file per sample, one variant per line (variable number of variants per sample). 
* column 1: variant key `chrom:pos:ref:alt` (`pos` is a 9-digit number with 0-prefix). 
* column 2: variant features (in JSON format)

```javascript
{
	"anno": {
		"af_1000g_all": float (0-1),
		"af_1000g_eas": float (0-1),
		"af_exac_all": float (0-1 or -1),
		"af_exac_adj": float (0-1 or -1),
		"af_exac_eas": float (0-1 or -1),
		"in_dbsnp": int (0/1),
		"in_exac_target": int (0/1),
	},
	"predict": {
		"pred": int (0/1),
		"prob": float (0-1),
	}
	"read": {
		"baseq_baq": list (11 int, 0-),
		"baseq_cal": list (11 int, 0-42),
		"baseq_raw": list (11 int, 0-42),
		"read_nt": list (11 int, 0-3),
		"uniq_24": list (12 float, 0-1, 1=uniq, 12=36-24),
		"map_pos": int (genome position),
		"mapq": int (0-42),
		"nsnp": int (1-2),
		"offset": int (0-31),
		"read_id": str ("tile_id:cluster_x:cluster_y"),
		"rev_comp": int (0/1), 
		"sample_id": str, 
	}
}
```

### allele frequency data

* data location: `data[234]/pacs/prov_af/[PROV]/af.b_[BID].sb_[SBID].out.gz` (Bin ID `BID` = `1..5733`, sub-bin ID `SBID` = `1..9`, `PROV` is the name of 34 provinces, including `All`, `Unknown` and `Korea`). 
* data size: 15.6TB.
* number of files: 51597 files per province (5733 bins, 9 sub-bins per bin). The number of variants per sub-bin is 500, except for the last bin (5733).  
* content: Bayesian allele frequency estimation result, one file per sub-bin, one variant per line. 
* column 1: variant key `chrom:pos:ref:alt` (`pos` is a 9-digit number with 0-prefix). 
* column 2: variant features (in JSON format)

```javascript
{
	"af_pdf": {
		"ci_mean_90%": list (2 int, 90% CI lower and upper bound around MLE),
		"ci_mean_95%": list (2 int, 95% CI lower and upper bound around MLE),
		"ci_mean_99%": list (2 int, 99% CI lower and upper bound around MLE),
		"ci_mean_99.5%": list (2 int, 99.5% CI lower and upper bound around MLE),
		"ci_mean_99.9%": list (2 int, 99.9% CI lower and upper bound around MLE),
		"ci_mean_99.99%": list (2 int, 99.99% CI lower and upper bound around MLE),
		"ci_mean_99.999%": list (2 int, 99.999% CI lower and upper bound around MLE),

		"ci_mode_90%": list (2 int, 90% CI lower and upper bound around mode),
		"ci_mode_95%": list (2 int, 95% CI lower and upper bound around mode),
		"ci_mode_99%": list (2 int, 99% CI lower and upper bound around mode),
		"ci_mode_99.5%": list (2 int, 99.5% CI lower and upper bound around mode),
		"ci_mode_99.9%": list (2 int, 99.9% CI lower and upper bound around mode),
		"ci_mode_99.99%": list (2 int, 99.99% CI lower and upper bound around mode),
		"ci_mode_99.999%": list (2 int, 99.999% CI lower and upper bound around mode),

		"init": {
			"af": float (0-1, initial value of AF),
			"index": int (0-999, index in AF search list),
			"prob": float (0-1),
		},
		"mean": {
			"af": float (0-1, MLE),
			"index": int (0-999, index in AF search list),
			"prob": float (0-1),
		},
		"median": {
			"af": float (0-1, median),
			"index": int (0-999, index in AF search list),
			"prob": float (0-1),
		},
		"mode": {
			"af": float (0-1, mode),
			"index": int (0-999, index in AF search list),
			"prob": float (0-1),
		},
		"pdf": [          # list of 1000 AF values and probabilities
			[af, prob],   # [float, float] 
			...
		], 
		"range": {
			"max": float (0-1, max AF search value), 
			"min": float (0-1, min AF search value), 
			"npoint": int (1000, number of AF search values), 
		},
	},
	"var_union": {
		"[DP]": {       # DP: int (0,1,2,...; 0=all depth),
			"aa": int (0,1,..., NS with only Alt allele),
			"ra": int (0,1,..., NS with both alleles),
			"rr": int (0,1,..., NS with only Ref allele),
			"ns": int (1,..., NS at current depth),
			"af": float (0-1, raw Alt allele frequency),
		},
		...
	},
}
```
### bam index (bai)

* data location: `data1/pacs/idx_stat/batch.[BID]/[SID].idxstats.gz` (Batch ID `BID`=`00..21`, `SID` is sample ID). 
* data size: 25MB
* number of files: 214277 (1 per sample, 10000 per batch)
* content: number of reads mapped to chr1-22, X, Y, MT (generated using `samtools idxstats [SID].bam`)
* columns (tab-separated): `chrom`, `length`, `#reads`, `#reads unmapped`(0).






