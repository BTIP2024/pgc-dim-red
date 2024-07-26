#' VCF to CSV conversion
#' 
#' Extract genotypes from VCF files and convert to CSV
#' Input can be a VCF file or excel file
#' @examples
#' input_file <- VCF_to_CSV("my_vcf_file.vcf");
#' @export

VCF_to_CSV <- function(input_file) {
  if(grepl('.vcf', input_file)){
    input <- vcfR::read.vcfR(input_file)
    input_geno <- extract.gt(input, return.alleles=TRUE)
    fixed_columns <- getFIX(input)
    input_final <- data.frame(fixed_columns, input_geno)
  } else {
    input_final <- readxl::read_excel(input_file)
  }
  
  readr::write_csv(input_final, file = "conversion_output.csv")
}