

count = 0
kegg <- readLines('~/Desktop/ko00001.keg')
kegg <- kegg[4:(length(kegg)-3)]



main=""
pathway_id=""
pathway_name=''
subpathway_id=''
subpathway_name=''
kegg_id=''
hgnc_id=''
note=''

dat <- data.frame(main="",pathway_id="",pathway_name='',subpathway_id='',subpathway_name='',kegg_id='',hgnc_id='',note='', stringsAsFactors = F)

lst <- list()


for (l in kegg){
  splitted <- unlist(strsplit(l, ';')) # unlist(strsplit(l, '\\ +'))
  first <- unlist(strsplit(splitted[1], '\\ +'))
  category <- first[1]
  if (unlist(strsplit(category, ''))[1] == 'A'){
    main = category
    pathway_id=""
    pathway_name=''
    subpathway_id=''
    subpathway_name=''
    kegg_id=''
    hgnc_id=''
    note=''
  }
  
  if (unlist(strsplit(category, ''))[1] == 'B'){
    nsplit <- unlist(strsplit(splitted, '\\ +'))
    pathway_id = nsplit[2]
    pathway_name = paste(nsplit[3:length(nsplit)],collapse=' ')
  }
  
  if (unlist(strsplit(category, ''))[1] == 'C'){
    nsplit <- unlist(strsplit(splitted, '\\ +'))
    subpathway_id = as.character(nsplit[2])
    subpathway_name = as.character(paste(nsplit[3:length(nsplit)],collapse=' '))
  }
  
  if (unlist(strsplit(category, ''))[1] == 'D'){
    kegg_id <- as.character(first[2])
    hgnc_id <- as.character(gsub(',', "",paste(first[3:length(first)], collapse = ';')))
    remain <- as.character(paste(splitted[2:length(splitted)], collapse = ' '))
  }
  
  #row=c(main,pathway_id,pathway_name ,subpathway_id ,subpathway_name ,kegg_id ,hgnc_id ,note )
  #dat <- rbind(dat, row)
  #print(row)
  count = count + 1
  print(round(count/length(kegg), 3))
  
  if (hgnc_id != ''){
    for (id in unlist(strsplit(hgnc_id, ';'))){
      if (is.null(lst[[id]])){
        lst[[id]] <- c(subpathway_name)
      }else{
        print('dis already here')
        lst[[id]] <- c(lst[[id]], c(subpathway_name))
      }
    }
  }

}






