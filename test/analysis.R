args <- commandArgs(trailingOnly = TRUE)
manual_data <- read.table(args[1],sep='\t',header=T)
my_data <- read.table(args[2],sep='\t',header=T)

manual_data$quant <- as.numeric(manual_data$quant)
manual_data$original_quant <- as.numeric(manual_data$original_quant)

my_data$annot <- my_data$quant
my_data$quant <- as.numeric(my_data$quant)


manual_data$key <- paste( manual_data$sequence, manual_data$glyco,sep=' ')
my_data$key <- paste( my_data$sequence, my_data$glyco,sep=' ')
my_singlets <- my_data[my_data$quant == 1e0-5 | my_data$quant == 1e05 | grepl("potential", my_data$annot ) | grepl("conflict",my_data$annot) ,]

merged <- merge(manual_data,my_singlets,by='key',all.x=T,all.y=T)
 
false_negatives <- subset(merged,!is.na(uniprot.x) & (quant.x == 0.00001 | quant.x == 100000) & is.na(uniprot.y))
positives <- subset(merged, ! is.na(uniprot.x) & quant.x == quant.y)
false_positives <- subset(merged, ! key %in% positives$key & ! is.na(uniprot.y) & !is.na(uniprot.x) & quant.y != quant.x & quant.y == original_quant & (original_quant == 1e05 | original_quant == 1e-05 ))
new_quants <- subset(merged, ! key %in% c( unique(positives$key), unique(false_positives$key) ) & ! is.na(quant.y) & original_quant == quant.y & (original_quant == 1e05 | original_quant == 1e-05))
invalid_negatives <- subset(merged, ! is.na(uniprot.x) & grepl("potential",annot) & quant.x < 100000 & quant.x > 0.00001 )
valid_negatives <- subset(merged, ! is.na(uniprot.x) & grepl("potential",annot) & (quant.x == 100000 | quant.x == 0.00001 ))
new_potentials <- subset(merged,is.na(quant.x) & grepl("potential",annot))

if (!is.na(args[3])) {
	write_excel(header=T,filename=args[3],
		false_negatives=false_negatives,
		false_positives=false_positives,
		positives=positives,
		new_singlets=new_quants,
		new_potentials=new_potentials,
		invalid_negatives=invalid_negatives,
		valid_negatives=valid_negatives)
}