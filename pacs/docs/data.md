This files documents the format of intermediate data in PACS.

### variant calling output 

* data location: `data2/pacs/var_json/batch_[BID]/[SID].var.json.gz`, (`BID` = `00..21`). 
* number of files: `214277`.
* content: Random Forest variant classification result, one file per sample, one variant per line. 
* column 1: variant key, format = `chrom:pos:ref:alt` (`pos` is a 9-digit number with 0-prefix). 
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
