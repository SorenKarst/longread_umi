# R-packages
require("Biostrings")
require("gridExtra")
require("tidyverse")
require("purrr")
require("ggrepel")
require("ggstance")
#require(ggsequence)

# Settings
options(scipen=999)
tt <- ttheme_default(core=list(fg_params=list(hjust=1, x = 0.95, fontsize = 7)),
                     colhead=list(fg_params=list(fontsize = 7)),
                     padding= unit(c(1.5, 1.5), "mm"))

# Format names
lu_name_format <- function(file, post_fix = ""){
  gsub(".*/|\\..*", "", file) %>%
    paste(., post_fix, sep = "")
}

# Calculate error rate from cigar
lu_cigar_err <- function(samdata){
  # Split cigar tags for each sequence
  cig_out <- str_extract_all(samdata$CIGAR, "\\d+|[MIDS]") %>%
    # Calculate MID and MD for each sequence
    {lapply(.,
            function(x){
              m <- matrix(x, ncol = 2, byrow = T)
              MID <- sum(
                as.integer(
                  m[,1][m[,2] != "S"]
                )
              )
              D <- sum(
                as.integer(
                  m[,1][m[,2] == "D"]
                )
              )
              I <- sum(
                as.integer(
                  m[,1][m[,2] == "I"]
                )
              )
              return(c(MID, I, D))
            }
    )} %>%
    # Convert til tibble
    do.call(rbind, .) %>%
    {tibble(MID = .[,1], I = .[,2], D = .[,3])}
  
  # Transform data
  samdata %>%
    transmute(umi = gsub(":.*", "", QNAME),
              target = RNAME,
              aln_len = cig_out$MID,
              error = as.integer(gsub("NM:i:", "", NM))/aln_len * 100,
              del = cig_out$D,
              ins = cig_out$I,
              mm = as.integer(gsub("NM:i:", "", NM)) - del - ins
    ) %>%
    ungroup %>%
    return()
}


lu_compile_qc <- function(
  umi_consensus,
  data_dir = "qc",
  reference = "zymo-ref-uniq_2019-03-15.fa",
  silva = "ssu_silva_db.fa", 
  read_data = "reads.fa",
  read_tax = "read_classification.txt",
  umi_bin_map = "umi_bin_map.txt",
  read_orientation = NULL,
  read_sam = lu_name_format(read_data,".sam"),
  ref_sam = lu_name_format(umi_consensus,".sam"),
  ref_ssu_sam = paste(lu_name_format(umi_consensus),
                      "_ssu_",
                      lu_name_format(reference),
                      ".sam",
                      sep =""),
  silva_ssu_sam = paste(lu_name_format(umi_consensus),
                        "_",
                        lu_name_format(silva, ".sam"),
                        sep =""),
  silva_tax = paste(lu_name_format(umi_consensus),
                    "_",
                    lu_name_format(silva, "_tax.txt"),
                    sep = ""),
  chimera = lu_name_format(umi_consensus, "_chimera.txt"),
  data_yield = "data_stats.txt",
  out_prefix = "",
  out_path = paste(
    "./",
    out_prefix,
    "qc.Rdata",
    sep = "")
){
  # Read data
  read_lengths <- Biostrings::fasta.seqlengths(
    paste(data_dir, "/", read_data, sep ="")
    )
  
  # Read taxonomy
  read_tax <- readr::read_delim(
    paste(data_dir, "/", read_tax, sep =""),
    delim = " ",
    col_names = F) %>%
    dplyr::transmute(read = X1,
              tax = gsub("_.*", "", X2)
              )
  
  # UMI read binning stats
  umi_bin_stats <- readr::read_delim(
    paste(data_dir, "/", umi_bin_map, sep =""),
    delim = " ",
    col_names = F) %>%
    dplyr::transmute(
      umi = gsub(";.*", "bins", X1),
      read = X2,
      umi_cluster_size = gsub(".*;", "", X1),
      umi_match_error = as.integer(X3)
      )
  
  umi_stats <- dplyr::left_join(umi_bin_stats, read_tax, by ="read") %>%
    group_by(umi, tax) %>%
    summarise(n = n(),
              umi_cluster_size = umi_cluster_size[1],
              umi_match_error = sum(umi_match_error)
              ) %>%
    group_by(umi) %>%
    summarise(umi_bin_size = sum(n),
              umi_cluster_size = as.integer(umi_cluster_size[1]),
              contamination = (sum(n) - max(n))/sum(n)*100,
              umi_match_error = sum(umi_match_error)/umi_bin_size)

  # UMI consensus length
  con_length <- Biostrings::fasta.seqlengths(
    paste(data_dir, "/", umi_consensus, sep ="")) %>%
    tibble::tibble(umi = gsub(";.*", "", names(.)),length = .)
  
  # Error
  sam_names = c("QNAME",
                "FLAG",
                "RNAME",
                "POS",
                "MAPQ",
                "CIGAR",
                "MNAME",
                "MPOS",
                "TLEN",
                "NM",
                "CS")
  
  ref_error <- readr::read_delim(
    file = paste(data_dir, "/", ref_sam, sep = ""),
    delim = "\t",
    col_names = sam_names) %>%
    lu_cigar_err() %>%
    dplyr::transmute(umi = gsub(";.*", "", umi),
                     ref_error = error,
                     ref_tax = target)
  
  if(!is.null(silva)){
    ref_ssu_error <- readr::read_delim(
      file = paste(data_dir, "/", ref_ssu_sam, sep = ""),
      delim = "\t",
      col_names = sam_names) %>%
      lu_cigar_err() %>%
      dplyr::transmute(umi = gsub(";.*", "", umi),
                       ref_ssu_error = error,
                       ref_ssu_tax = target)
    
    silva_error <- readr::read_delim(
      file = paste(data_dir, "/", silva_ssu_sam, sep = ""),
      delim = "\t",
      col_names = sam_names) %>%
      lu_cigar_err() %>%
      dplyr::transmute(umi = gsub(";.*", "", umi),
                       silva_ssu_error = error,
                       silva_ssu_id = target)
    
    silva_tax <- readr::read_delim(
      paste(data_dir, "/", silva_tax, sep =""),
      delim = " ",
      col_names = c("silva_ssu_id", "silva_tax")) %>%
      tidyr::separate(silva_tax,
                      c("silva_k", "silva_p", "silva_c", "silva_o", "silva_f", "silva_g", "silva_s"),
                      sep = ";", 
                      fill = "right")
    
    silva <- silva_error %>%
      left_join(silva_tax, by = "silva_ssu_id") 
  }
  
  # Raw read orientation ratio (+/-)
  if(!is.null(read_orientation)){  
    ror <- readr::read_delim(
      file = paste(data_dir, "/", read_orientation, sep = ""),
      delim = " ",
      col_names = c("umi", "+", "-", "?")) %>%
      dplyr::transmute(umi,
                       ror = `+`/`-`
                       )
  }
  
  # Chimera test
  chimera <- readr::read_delim(
    paste(data_dir, "/", chimera, sep =""),
    delim = "\t",
    col_names = F) %>%
    dplyr::transmute(umi = gsub(";.*", "", X1),
                     chimera = X3)    

  # Yield stats  
  yield_stats <- readr::read_delim(
    paste(data_dir, "/", data_yield, sep = ""),
    delim = ","
  )
  
  # Ref homopolymer statistics
  ref_hp <- Biostrings::readDNAStringSet(
    paste(data_dir, "/", reference, sep = "")
  ) %>%
    lu_ref_hp_pos()
  
  # Check map efficiency
  map_check <- left_join(con_length,
                         ref_error,
                         by = "umi") %>%
    summarise(Mapped = sum(!is.na(ref_error)),
              Unmapped = sum(is.na(ref_error)))
  
  # Read error
  read_error <- readr::read_delim(
    file = paste(data_dir, "/", read_sam, sep = ""),
    delim = "\t",
    col_names = sam_names) %>%
    lu_cigar_err() %>%
    dplyr::transmute(umi = gsub(";.*", "", umi),
                     ref_error = error,
                     ref_tax = target)
  
  # Combine data
  qc <-  umi_stats %>%
    left_join(con_length, by = "umi") %>%
    left_join(ref_error, by = "umi") %>%
    {if(!is.null(silva))left_join(., ref_ssu_error, by = "umi") else .} %>%
    left_join(chimera, by = "umi") %>%
    {if(!is.null(silva))left_join(., silva, by = "umi") else .} %>%
    {if(!is.null(read_orientation))left_join(., ror, by = "umi") else .}
    
  
  # Prepare output
  out <- c(
    qc = paste(out_prefix, "qc", sep = ""),
    read_lengths = paste(out_prefix, "read_lengths", sep = ""),
    yield_stats = paste(out_prefix, "yield_stats", sep = ""),
    map_check = paste(out_prefix, "map_check", sep = ""),
    ref_hp = paste(out_prefix, "ref_hp", sep = ""),
    read_error = paste(out_prefix, "read_error", sep = "")
  )
  
  # Assign with prefix for export
  assign(out["qc"], qc)
  assign(out["read_lengths"], read_lengths)
  assign(out["yield_stats"], yield_stats)
  assign(out["map_check"], map_check)
  assign(out["ref_hp"], ref_hp)
  assign(out["read_error"], read_error)

  
  save(
    list = out,
    file = out_path
  )
}

# Plot histogram
lu_plot_hist <- function(value, bin_width, fill = NULL){
  ggplot(tibble(value), aes(value, fill = fill)) +
    geom_histogram(binwidth = bin_width)
}

# Plot yield stats
lu_yield_summary <- function(
  yield_stats,
  qc,
  title){
  # Calculate reads count and bp
  process_stats <- tibble(
    data_type = c(
      "UMI clusters",
      "Binned",
      "Assembled"),
    read_count = c(
      sum(as.integer(qc$umi_cluster_size), na.rm = T),
      sum(qc$umi_bin_size, na.rm = T),
      filter(qc, !is.na(length)) %>% {sum(.$umi_bin_size, na.rm = T)}
      ),
    bp_average = mean(qc$length, na.rm =T)
    ) %>%
    mutate(bp_total = read_count * bp_average,
           bins_n = c("NA", nrow(qc), filter(qc, !is.na(length)) %>% nrow()))

  # Combine stats
  stats <- bind_rows(yield_stats, process_stats) %>%
    mutate(`%_bp` = round(bp_total/bp_total[1]*100))
  colnames(stats) <- c(
    "Processing step",
    "Read count",
    "Total bp",
    "Average bp",
    "UMI bin\ncount",
    "% input bp"
  )
  # Plot stats table
  grid.arrange(tableGrob(stats, rows= NULL, theme = tt), top = title)
}

# Plot relative abundance skew as function of operons
lu_binsize_skew_plot<- function(qc, top_n = 10){
  # Calculate skew
  skew <- qc %>%
    drop_na() %>%
    group_by(ref_tax) %>%
    summarise(median = median(umi_bin_size),
              sd = sd(umi_bin_size)) %>%
    mutate(sd = if_else((is.na(sd) | sd == 0), median, sd)) %>%
    rowwise() %>%
    mutate(skew = as.character(cut(median,
                                   breaks = c(0, sd*2, sd*3, Inf),
                                   labels = c("danger", "warning", "good"),
                                   include.lowest = T
    )),
    skew_ratio = median/sd) %>%
    as_tibble() %>%
    arrange(skew_ratio) %>%
    mutate(ref_tax_grp = gsub("_.*", "", ref_tax))
  
  # Plot skew catergory counts 
  sp1 <- skew %>%
    group_by(ref_tax_grp, skew) %>%
    summarise(n = n())
  
  # Plot top skew densities
  skew_top <- skew %>%
    filter(skew %in% c("danger", "warning")) %>%
    top_n(- top_n, wt = skew_ratio)

  sp2 <- ggplot(filter(qc, ref_tax %in% skew_top$ref_tax),
                aes(umi_bin_size,
                    stat(count),
                    color = gsub("_.*", "", ref_tax),
                    group = ref_tax)) +
    geom_density() +
    scale_color_discrete(name = "Operon") +
    facet_grid(rows = vars(gsub("_.*", "", ref_tax)))
  
  # Generate plot
  grid.arrange(ggplotGrob(sp2),
               tableGrob(sp1, rows= NULL, theme = tt),
               nrow = 2)
}

# Plot chimera overview table
lu_chimera_summary <- function(qc, min_bin_size = 10){
  # Summarise chimera stats
  chimera <- qc %>%
    filter(umi_bin_size >= min_bin_size) %>%
    summarise(`Chimera-` = sum(chimera == "N", na.rm = T),
              `Chimera+` = sum(chimera == "Y", na.rm = T),
              `Chimera?` = sum(chimera == "?", na.rm = T))
  
  # Plot chimera table
  grid.arrange(tableGrob(chimera, rows= NULL, theme = tt))
}

# Plot error scatter plot
lu_error_plot <- function(qc,
                          x = "umi_bin_size",
                          y = "ref_error",
                          min_bin_size = 0,
                          point_size = 1,
                          flags = F){
  # Filter sequences
  qcf <- filter(qc,
                umi_bin_size >= min_bin_size,
                !is.na(length))
  
  # Build plot
  pe <- ggplot(qcf,
               aes_string(x = x,
                          y = y)
  ) +
    geom_point(size = point_size, color = alpha("black", 0.1)) +
    theme_classic() +
    theme(legend.position = "none",
          text = element_text(size = 8),
          axis.text = element_text(size = 6))
  if (flags){
    pe <- pe +
      geom_point(size = point_size, aes(color = contamination)) +
      geom_point(size = 3,
                 shape = 21, 
                 fill = NA,
                 color = if_else(qcf$chimera == "Y", "Blue", "NA")) +
      scale_color_gradient(low = NA, high = alpha("red", 1))
  }
  
  # Plot
  pe
}


# Plot error table
lu_error_plot_tbl <- function(qc,
                              x = "umi_bin_size",
                              y = "ref_error",
                              breaks = c(seq(0,100,5), Inf),
                              min_bin_size = 0,
                              digits = 4){
  # Filter sequences
  qcf <- filter(qc,
                umi_bin_size >= min_bin_size,
                !is.na(length))
  
  # Calculate error pr. bin
  qcf_tbl <- data.frame(x = qcf[,x] %>% unlist(),
                        y = qcf[,y] %>% unlist()) %>%
    mutate(grp = cut(x, breaks)) %>%
    group_by(grp) %>%
    summarise(y = round(mean(y), digits), n = n())
  colnames(qcf_tbl) <- c(paste(x, "grp", sep ="_"), y, "n")
  
  # Plot
  grid.arrange(tableGrob(qcf_tbl, rows= NULL, theme = tt))
}

# Plot variant error table
lu_variant_error_tbl <- function(
  profile,
  title,
  min_bin_size = 0,
  digits = 4
){
  # Format data
  pf <- profile %>%
    mutate(
      rname = gsub(";.*", "", rname),
      umi_bin_size = as.integer(gsub(".*=", "", qname)),
      qname = gsub(";.*", "", qname)
    ) %>%
    filter(
      umi_bin_size >= min_bin_size,
      !is.na(ref_alnlen)
    )
  
  # Summarise errors
  pf_tbl <- profile %>%
    mutate(
      rname = gsub(";.*", "", rname),
      cluster_size = as.integer(gsub(".*=", "", qname)),
      qname = gsub(";.*", "", qname)
    ) %>%
    group_by(qname, ref_alnlen, rname, category, cluster_size) %>%
    summarise(error = sum(count)) %>%
    spread(key = category, value = error, fill = 0) %>%
    ungroup() %>%
    # Necessary if some error types are not present.
    bind_rows(
      tibble(
        qname = as.integer(),
        ref_alnlen = as.integer(),
        rname = as.integer(),
        cluster_size = as.integer(),
        mm = as.numeric(),
        del = as.numeric(),
        ins = as.numeric(),
        perfect = as.numeric()
      )
    ) %>%
    mutate_at(vars(mm, ins, perfect, del), ~if_else(is.na(.), as.numeric(0), .)) %>%
    mutate(error = round((del+ins+mm)/ref_alnlen*100, digits)) %>%
    arrange(gsub("_.*", "", rname), as.integer(gsub(".*_|;.*", "", rname)), error) %>%
    select(rname, qname, del, ins, mm, error, cluster_size) %>%
    group_by(rname) %>%
    mutate(
      type = if_else(
        as.integer(cluster_size) == max(as.integer(cluster_size)),
        "Best",
        "Spurious"
      )
    )
}

# Plot variant error scatter plot
lu_variant_error_plot <- function(
  var_table,
  title = NULL
){
  # Format variant data
  vtf <- var_table %>% 
    ungroup() %>%
    mutate(
      species = gsub("_.*", "", rname),
      operon_n = gsub("^[a-zA-Z]*_|;.*|_.*", "", rname) %>%
        as.integer()
    ) %>%
    arrange(operon_n) %>%
    mutate(
      operon = as_factor(gsub("^[a-zA-Z]*_|;.*", "", rname))
    )

  #Plot
  ggplot(
    vtf,
    aes(
      x=operon,
      y = del+ins+mm,
      color = species,
      size = cluster_size
    )
  ) +
    geom_point() +
    #scale_y_continuous(limits = c(0, 10)) +
    scale_size(range = c(0.1, 5)) +
    facet_grid(~species, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = 'off') +
    theme_classic() +
    labs(
      x = "Operon (#)",
      y = "Errors (n)",
      title = title
    ) +
    theme_classic() +
    theme(strip.background = element_blank(),
          legend.position = "none",
          text = element_text(size = 6),
          axis.text = element_text(size = 4),
          axis.text.x = element_text(angle = 45),
          axis.ticks = element_line(size = 0.25, color = "darkgrey"),
          axis.ticks.length = unit(2, "pt"),
          axis.line = element_line(size = 0.25)) +
    annotate(
      "segment",
      x=Inf,
      xend=-Inf,
      y=Inf,
      yend=Inf,
      color="black",
      lwd=1
    )  
}


# Plot contamination overview statistics 
lu_contamination_plot <- function(qc, min_bin_size = 10){
  # Summarise contamaintion
  contamination <- qc %>%
    filter(umi_bin_size >= min_bin_size, !is.na(length)) %>%
    mutate(contamination_grp = cut(contamination, seq(0, 100, 5), include.lowest = T)) %>%
    group_by(contamination_grp) %>%
    summarise(mean_error = mean(ref_error, na.rm = T), n = n())
  
  # Plot contamiantion table
  grid.arrange(tableGrob(contamination, rows= NULL, theme = tt))
}

# Model error and flag artefacts
lu_artefact_plot <- function(qc,
                             breaks = c(seq(0, 50, 3), Inf),
                             iqr_prob = 0.9,
                             offset = 0.25,
                             sd_x = 7,
                             iqr = NULL,
                             out_type= "plot",
                             point_size = 1){
  # Filter data
  qcf <- filter(qc,
                !is.na(length))
  
  # Define UMI bin size intervals
  qc_ubsi <- qcf %>%
    mutate(ubsi = cut(umi_bin_size, breaks, include.lowest = T),
           ref_error_n = ceiling(ref_error * length/100))
  
  # Calculate quantiles ranges and cutoffs
  if (is.null(iqr)){
    iqr <- qc_ubsi %>%
      group_by(ubsi) %>%
      summarise(iqr_n = quantile(ref_error_n,
                               probs = iqr_prob,
                               na.rm = T,
                               type = 8),
                median_n = median(ref_error_n),
                cutoff_n = median_n + iqr_n,
                mean_ubs = mean(umi_bin_size),
                n = n()
      )    
  }
  if (out_type == "iqr"){
    return(iqr)
  }
  
  # Trim data
  qc_ubsi_f <- left_join(qc_ubsi, iqr, by = "ubsi") %>%
    filter(ref_error_n <= cutoff_n)
  
  # Fit model to trimmed data
  err_mod <- mgcv::gam(
    formula = ref_error_n ~ s(umi_bin_size, bs = "cs"),
    family = "poisson",
    data = qc_ubsi_f)
  
  if (out_type == "model"){
    return(model)
  }
  
  # Calculate cutoff
  qc_ubsi_f <- qc_ubsi_f %>%
    mutate(pred = mgcv::predict.gam(object = err_mod,
                                    type="response")/length*100
    )
  
  qc_co <- qc_ubsi_f %>%
    mutate(pred_err = pred - ref_error) %>%
    group_by(ubsi) %>%
    summarise(mean_ubs = mean(mean_ubs),
              mean_pred = mean(pred),
              sd = sd(pred_err),
              cutoff = mean_pred + offset + sd*sd_x) %>%
    add_row(ubsi = "inf",
            mean_ubs = max(qcf$umi_bin_size),
            mean_pred = .$mean_pred[nrow(.)],
            sd = .$sd[nrow(.)],
            cutoff = .$cutoff[nrow(.)])
  flag <- left_join(qc_ubsi, qc_co, by = "ubsi") %>%
    mutate(aflag = if_else(ref_error > cutoff, "lightblue","NA"))
  
  # Plot
  lu_error_plot(flag, point_size = point_size) +
    geom_point(color = flag$aflag, size = point_size) +
    geom_line(data = qc_co, aes(mean_ubs, mean_pred), color = "red") +
    geom_line(data = qc_co, aes(mean_ubs, cutoff), color = "blue")
}

# Determine homopolymer positions in references
lu_ref_hp_pos <- function(refs){
  lapply(
    refs,
    function(seq){
      seq <- as.character(seq)
      homopolymer = str_extract_all(seq, "A{3,}|T{3,}|C{3,}|G{3,}") %>%
        unlist()
      startpos = gregexpr("A{3,}|T{3,}|C{3,}|G{3,}", seq)[[1]]
      hlen = paste("H",nchar(homopolymer), sep = "")
      hbase = substring(homopolymer, 1, 1)
      stoppos = as.integer(startpos + nchar(homopolymer)-1)
      
      # Output
      return(tibble(homopolymer, hbase, startpos, stoppos, hlen))
    }
  )
}

# Convert CS tag from sam file
lu_parse_sam_cs <- function(input) {
  
  # Align length function
  fun_alnlen <- function(cigar_entry){
    # Split cigar tags for each sequence
    str_extract_all(cigar_entry, "\\d+|[MIDS]") %>%
      # Calculate MID and MD for each sequence
      {lapply(.,
              function(x){
                m <- matrix(x, ncol = 2, byrow = T)
                MID <- sum(
                  as.integer(
                    m[,1][m[,2] != "S"]
                  )
                )
                MD <- sum(
                  as.integer(
                    m[,1][m[,2] %in% c("M", "D")]
                  )
                )
                return(c(MID, MD))
              }
      )} %>%
      # Convert til tibble
      do.call(rbind, .) %>%
      {tibble(MID = .[,1], MD = .[,2])}
  }
  
  # Calculate align lengths
  inputf <- input %>%
    bind_cols(fun_alnlen(.$CIGAR))

  # Process error blocks
  err_blocks <-
    # Split CS tags into components for each sequence
    str_extract_all(inputf$CS, ":\\d+|[-+*][actg]{1,}") %>%
    # Map deletions to individual positions
    lapply(
      function(x){
        l <- as.list(x)
        if_else(
          grepl("^-", l),
          str_sub(l, 2) %>%
            str_split("") %>%
            lapply(., function(x)gsub("^", "-",x)),
          l
        ) %>%
          unlist() %>%
          list(value = .)
      } 
    ) %>%
    # Convert data to tibble
    setNames(., inputf$QNAME) %>%
    bind_rows(.id = "QNAME") %>%
    # Profiling
    mutate(
      # Error type
      type = if_else(grepl("\\:", value), NA_character_, value),
      category = str_extract(value, "[:*\\-+]"),
      category = case_when(
        category == ":" ~ "perfect",
        category == "*" ~ "mm",
        category == "-" ~ "del",
        category == "+" ~ "ins",
      ),
      # Error lengths in query and ref space
      count = case_when(
        category == "perfect" ~ sub("\\:| ", "", value),
        category == "mm" ~ as.character(str_count(sub("\\*", "", value))/2),
        category == "del" ~ as.character(str_count(sub("\\-", "", value))),
        category == "ins" ~ as.character(str_count(sub("\\+", "", value))),
      ),
      count = as.integer(count),
      ref_count = case_when(
        category == "ins" ~ as.integer(0),
        TRUE ~ count
      )
    ) %>%
    # Calculate positions relative to reference
    group_by(QNAME) %>%
    mutate(
      ref_stop = cumsum(ref_count),
      ref_start = if_else(
        category == "ins",
        as.integer(ref_stop - ref_count),
        as.integer(ref_stop - ref_count + 1)
      )
    )
  
  # Merge sam data with error blocks
  sam_error <- right_join(
    inputf[c("QNAME", "RNAME", "FLAG", "POS", "MID", "MD")],
    err_blocks[c("QNAME", "category", "type", "ref_start", "count")],
    by = "QNAME"
  ) %>%
    transmute(
      qname = QNAME,
      rname = RNAME,
      flag = FLAG,
      alnlen = MID,
      ref_alnstart = POS,
      ref_alnlen = MD,
      category,
      type,
      ref_pos = ref_start + ref_alnstart - 1,
      count
    )
  
  # Convert to MA style output
  sam_error <- sam_error %>%
    mutate(
      count = if_else(category == "perfect", as.integer(0), count),
    ) %>%
    group_by(qname) %>%
    filter(!(category == "perfect" & n() > 1))
  
  # Output
  return(sam_error)
}

# Generate error profiles
lu_error_profile <- function(refs, sam){
  # Import data
  sam_names = c("QNAME",
                "FLAG",
                "RNAME",
                "POS",
                "MAPQ",
                "CIGAR",
                "MNAME",
                "MPOS",
                "TLEN",
                "NM",
                "CS")
  
  # Sam data
  if (is.character(sam)){
    sam <- read_delim(file = sam,
                      delim = "\t",
                      col_names = sam_names)
    
    # Check CS tag integrity
    sam_check <- !grepl("cs:Z:", sam$CS)
    if (sum(sam_check > 0)){
      message(
        paste(
          "Warning: ",
          sum(sam_check),
          " entries did not contain `cs:Z:` in CS column:",
          sep = ""
        )
      )
      print(sam[sam_check, ])
      message(
        "The entries have been removed."
      )
      sam <- sam[!sam_check, ]
    }
  }
  
  # Ref sequences
  if (is.character(refs)){
    refs <- readDNAStringSet(filepath = refs)
  }
  
  # Detect homopolymer regions
  homo <- lu_ref_hp_pos(refs)
  homo_ranges <- homo %>%
    lapply(
      function(ref){
        #Name ranges
        ref <- ref %>%
          mutate(
            range_name = paste(startpos, stoppos, sep = "-")
          )
        
        # Expand ranges
        ref_exp <- ref %>%
          select(range_name, startpos, stoppos) %>%
          apply(
            1,
            function(r){
              seq.int(r["startpos"], r["stoppos"])
            }
          ) %>%
          setNames(paste(ref$range_name, "_", sep = ""))%>%
          {tibble(
            range_name = names(unlist(.)),
            ref_pos = unlist(.)
          )} %>%
          mutate(range_name = gsub("_.*", "", range_name)) %>%
          right_join(ref, ., by = "range_name")
      }
    ) %>%
    bind_rows(.id = "rname") %>%
    select(rname, homopolymer, hpos = startpos, ref_pos)
  
  # Profile errors
  sam_errors <- lu_parse_sam_cs(sam)  

  # Link errors with homopolymer regions
  err_profile <- sam_errors %>% 
    ungroup() %>%
    left_join(
      homo_ranges,
      by = c("rname", "ref_pos")
    )

  # Output
  return(err_profile)
} 


# Plot error positions on sequences
lu_errorpos_plot <- function(profile,
                               qc,
                               species,
                               labels = NULL,
                               flag = NULL
                             ){
  # Format profile
  pf <- profile %>%
    separate(rname, c("pspecies", "poperon"), sep = "_", extra = "drop", remove = F) %>%
    mutate(qname = gsub(";.*", "", qname),
           pstart = 1,
           pend = ref_alnlen,
           perr = ref_pos - ref_alnstart + 1,
           poperon = gsub("_.*|;.*", "", poperon),
           poperon_dodge <- as_factor(poperon)) %>%
    filter(perr >= pstart & perr <= pend,
           grepl(species, pspecies),
           category != "perfect",
           qname %in% qc$umi)
  
  if (nrow(pf) > 0){
    # Estimate operon lengths
    pref <- pf %>%
      distinct(qname, .keep_all = T) %>%
      transmute(qname,
                pspecies,
                poperon,
                pstart,
                pend) %>%
      gather(key = pos,
             value = perr,
             pstart:pend)
    
    # Baseplot
    p <- ggplot(pf,
                aes(x = perr, y = qname)) +
      geom_line(data = pref, color = "gray", linetype = "dashed") +
      geom_point(data = pref, color = "gray") +
      geom_point(aes(color = category), position=ggstance::position_dodgev(height=0.3)) +
      facet_grid(rows = vars(poperon),
                 cols = vars(pspecies),
                 scales = "free",
                 space = "free_y",
                 drop = T) +
      theme_classic() +
      theme(text = element_text(size = 8),
            axis.text = element_text(size = 6))
    
    # Add labels
    if (!is.null(labels)){
      p <- p + 
        geom_text_repel(aes(label = perr))
      
    }
    
    # Add flags
    if (!is.null(flag)){
      # Format flags
      qcf <- qc  %>%
        separate(ref_tax,
                 c("pspecies", "poperon"),
                 sep = "_",
                 extra = "drop",
                 remove = F) %>%
        mutate(qname = umi,
               perr = 0,
               poperon = gsub("_.*|;.*", "", poperon)) %>%
        filter(grepl(species, pspecies))
      
      # Add flags to plot
      p <- p + 
        geom_text(data = qcf,
                  aes_string(x = -500,
                             y = "qname",
                             label = flag),
                  size = 2,
                  hjust = 0.3,
                  vjust = 0.1)
    }
    
    # Plot
    p
    
  } else {
    print(paste("No sequences of ",
                species,
                " with errors in current dataset.", sep = ""))
  }
}

# Function for calculating homopolymer fraction in alignment length
lu_hp_len <- function(qname, 
                       rname,
                       aln_start,
                       aln_stop,
                       ref_hp){
  
  # Calculate hp fraction in alignment pr. sequence
  purrr::pmap_dfr(
    list(qname, rname, aln_start, aln_stop, ref_hp = list(ref_hp)),
    function(qname, rname, aln_start, aln_stop, ref_hp){
      ref_hp[[rname]] %>%
        filter(
          startpos >= aln_start,
          stoppos <= aln_stop) %>%
        group_by(homopolymer) %>%
        summarise(
          hl_n = n(),
          hl_total = sum(stoppos - startpos + 1) 
        ) %>%
        mutate(qname)
    }
  )
}

# Error type line as a function of UMI bin size plot
lu_errortype_plot <- function(
  profile,
  ref_hp,
  break_size = 5,
  min_size = 3,
  legend_pos = "none"
){
  # Extract overall alignment data from each sequence
  aln_meta <- profile %>%
    distinct(
      qname,
      rname,
      ref_alnstart,
      ref_alnlen
      )
  
  # Calculate hp+ and hp- lengths for each sequence
  hp_len <- aln_meta %>%
    {lu_hp_len(
      .$qname,
      .$rname,
      .$ref_alnstart,
      .$ref_alnlen + .$ref_alnstart - 1,
      ref_hp
    )} %>%
    group_by(qname) %>%
    summarise(`hp+` = sum(hl_total)) %>%
    left_join(aln_meta[c("qname", "ref_alnlen")], by = "qname") %>%
    transmute(qname, `hp+`, `hp-` = ref_alnlen - `hp+`) %>%
    gather(
      key = "region_type",
      value = "region_length",
      `hp+`,
      `hp-`
    )

  
  # Summarise based on error type
  error_type <- profile %>%
    mutate(region_type = if_else(is.na(homopolymer), "hp-", "hp+")) %>%
    group_by(qname, category, region_type ) %>%
    summarise(count = sum(count)) %>%
    filter(category != "perfect") %>%
    left_join(
      tibble(
        qname = rep(aln_meta$qname, each = 2 * 3),
        region_type = rep(rep(c("hp-", "hp+"), each = 3), times = length(aln_meta$qname)),
        category = rep(c("mm", "ins", "del"), times = 2 * length(aln_meta$qname))
      ), 
      ., 
      by = c("qname", "region_type", "category")
    ) %>%
    ungroup() %>%
    mutate(count = if_else(is.na(count), as.integer(0), count)) %>%
    left_join(hp_len, by = c("qname", "region_type")) %>%
    mutate(
      umi_bin_size = as.integer(gsub(".*ubs=", "", qname)),
      qname = gsub(";.*", "", qname),
      ubs_grp = {
        breaks= c(seq(0,100,break_size)-1, Inf)
        lab = paste(
          c(min_size,seq(ceiling(min_size/break_size)*break_size, 100 - break_size, break_size), ""),
          c(paste("-", seq(ceiling(min_size/break_size)*break_size, 100,break_size)-1, sep = ""),"100+"),
          sep = ""
        )
        cut(umi_bin_size, breaks, lab)
      }
    ) %>%
    group_by(region_type, category, ubs_grp) %>%
    summarise(mean_err = sum(count)/sum(region_length)*100)
  
  # Plot
  ggplot(error_type, 
         aes(
           x = ubs_grp,
           y = mean_err,
           group = paste(region_type, category), 
           color = category, linetype = region_type
         )
  ) +
    geom_line(size = 0.5) +
    theme_classic() +
    labs(colour = "Error type",
         x = "UMI consensus read coverage (bins)",
         y = "Mean error-rate (%)",
         linetype = "Region type") +
    theme(
      text = element_text(size = 6),
      axis.text = element_text(size = 6),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
      legend.spacing.y = unit(0, 'lines'),
      legend.position = legend_pos,
      legend.title = element_text(size = 4),
      legend.text = element_text( size = 4),
      legend.key = element_rect(size = 1, linetype = 0),
      legend.key.size = unit(0.5, 'lines')
    ) +
    guides(
      shape = guide_legend(override.aes = list(size = 0.1)),
      color = guide_legend(override.aes = list(size = 0.1))
    )
}

# Plot errortype as function of UMI bin size table
lu_errortype_plot_tbl <- function(
                              profile,
                              break_size = 5,
                              digits = 3,
                              title){
  
  # Summarise based on error type
  errortype_tbl <- profile %>%
    # Assign bin size groups
    mutate(
      umi_bin_size = as.integer(gsub(".*ubs=", "", qname)),
      qname = gsub(";.*", "", qname),
      ubs_grp = {
        breaks= c(seq(0,100,break_size)-1, Inf)
        lab = paste(
          c("2",seq(break_size,100 - break_size, break_size), ""),
          c(paste("-", seq(break_size,100,break_size)-1, sep = ""),"100+"),
          sep = ""
        )
        cut(umi_bin_size, breaks, lab)
      }
    ) %>%
    # Sum errors pr. sequence
    group_by(qname, ref_alnlen, category, ubs_grp) %>%
    summarise(
      count = sum(count)
    ) %>%
    # Convert erros to wide form
    spread(
      key = category,
      value = count
    ) %>%
    # Calculate rates pr UMI bin size group
    group_by(ubs_grp) %>%
    summarise(
      all_err = sum(del, ins, mm, na.rm = T)/sum(ref_alnlen)*100,
      mm_err = sum(mm, na.rm = T)/sum(ref_alnlen)*100,
      ins_err = sum(ins, na.rm = T)/sum(ref_alnlen)*100,
      del_err = sum(del, na.rm = T)/sum(ref_alnlen)*100,
      seq_n = n()
    ) %>%
    # Format decimals
    mutate_at(
      vars(all_err:del_err),
      ~round(., digits)
    )
    colnames(errortype_tbl) <- c(
      "UMI bin size",
      "Total\nError rate (%)",
      "Mismatch\nError rate (%)",
      "Insert\nError rate (%)",
      "Deletion\nError rate (%)",
      "Sequence\n count (n)"
      )
  
  # Plot
  grid.arrange(tableGrob(errortype_tbl, rows= NULL, theme = tt), top = title)
}

# Error type summary statistics
lu_errortype_summary <- function(profile,
                                 title,
                                 digits=3){
  
  # Extract overall alignment data from each sequence
  aln_meta <- profile %>%
    distinct(
      qname,
      rname,
      ref_alnstart,
      ref_alnlen
    )
  
  # Calculate hp+ and hp- lengths for each sequence
  hp_len <- aln_meta %>%
    {lu_hp_len(
      .$qname,
      .$rname,
      .$ref_alnstart,
      .$ref_alnlen + .$ref_alnstart - 1,
      ref_hp
    )} %>%
    group_by(qname) %>%
    summarise(`hp+` = sum(hl_total)) %>%
    left_join(aln_meta[c("qname", "ref_alnlen")], by = "qname") %>%
    transmute(qname, `hp+`, `hp-` = ref_alnlen - `hp+`) %>%
    gather(
      key = "region_type",
      value = "region_length",
      `hp+`,
      `hp-`)
    
  
  # Summarise based on error type
  error_type <- profile %>%
    mutate(region_type = if_else(is.na(homopolymer), "hp-", "hp+")) %>%
    group_by(qname, category, region_type) %>%
    summarise(count = sum(count)) %>%
    # Convert erros to wide form
    spread(
      key = category,
      value = count
    ) %>%
    left_join(hp_len,., by = c("qname", "region_type"))
  
  errortype_tbl1 <- error_type %>%
    ungroup() %>%
    summarise(
      `Region type` = "all",
      `Mismatch\nError rate (%)` = sum(mm, na.rm = T)/sum(region_length)*100,
      `Insert\nError rate (%)` = sum(ins, na.rm = T)/sum(region_length)*100,
      `Deletion\nError rate (%)` = sum(del, na.rm = T)/sum(region_length)*100,
      `Total\nError rate (%)` = sum(mm, ins, del, na.rm = T)/sum(region_length)*100
    )
  
  errortype_tbl2 <- error_type %>%
    # Calculate rates pr UMI bin size group
    group_by(region_type) %>%
    summarise(
      `Mismatch\nError rate (%)` = sum(mm, na.rm = T)/sum(region_length)*100,
      `Insert\nError rate (%)` = sum(ins, na.rm = T)/sum(region_length)*100,
      `Deletion\nError rate (%)` = sum(del, na.rm = T)/sum(region_length)*100,
      `Total\nError rate (%)` = sum(mm, ins, del, na.rm = T)/sum(region_length)*100
    ) %>%
    dplyr::rename(`Region type` = region_type) %>%
    bind_rows(
      errortype_tbl1
    ) %>%
    # Format decimals
    mutate_at(
      vars(`Mismatch\nError rate (%)`:`Total\nError rate (%)`),
      ~round(., digits)
    )

  # Plot
  grid.arrange(tableGrob(errortype_tbl2, rows= NULL, theme = tt), top = title)   
}

# Error type as function of homopolymer
lu_hperror_summary <- function(
  profile,
  ref_hp,
  title
){
  
  # Extract overall alignment data from each sequence
  aln_meta <- profile %>%
    distinct(
      qname,
      rname,
      ref_alnstart,
      ref_alnlen
    )
  
  # Calculate hp+ and hp- lengths for each sequence
  hp_len <- aln_meta %>%
    {lu_hp_len(
      .$qname,
      .$rname,
      .$ref_alnstart,
      .$ref_alnlen + .$ref_alnstart - 1,
      ref_hp
    )} %>%
    group_by(homopolymer) %>%
    summarise(
      hp_len = sum(hl_total),
      hp_n = sum(hl_n)
    )
  
  # Summarise based on homopolymer
  hp_error <- profile %>%
    filter(!is.na(homopolymer)) %>%
    group_by(homopolymer) %>%
    summarise(error_n = sum(count)) %>%
    left_join(hp_len, ., by = "homopolymer") %>%
    mutate(
      error_rate = error_n/hp_len*100,
      hp_type = substr(homopolymer, 1, 1),
      hp_len = str_count(homopolymer)
    ) %>% mutate_at(
      vars(error_n:error_rate),
      ~replace_na(., 0)
    )
  
  # Error rate table
  error_tbl <- hp_error %>%
    select(hp_len, hp_type, error_rate) %>%
    spread(
      key = hp_type,
      value = error_rate, 
      fill = NA
    ) %>%
    mutate_at(
      vars(`A`:`T`),
      ~round(., 3)
    ) %>%
    mutate_all(
      as.character
    )
  
  # Count table
  n_tbl <- hp_len %>%
    transmute(
      hp_len = str_length(homopolymer),
      hp_type = str_sub(homopolymer, 1, 1),
      hp_n
    ) %>%
    select(hp_len, hp_type, hp_n) %>%
    spread(
      key = hp_type,
      value = hp_n, 
      fill = 0
    ) %>%
    mutate_all(
      as.character
    )
  
  
  # Merge tables
  tbl <- bind_rows(
    error_tbl,
    tibble(
      hp_len = "",
      A = "",
      `T` = "",
      C = "",
      G = ""
    ),
    n_tbl
  )
  
  # Plot
  grid.arrange(tableGrob(tbl, rows= NULL, theme = tt), top = title)   
}

# Plot error frequency as function of reference position
lu_ref_error_frequency <- function(
  profile,
  species = "",
  lower_cut_off = 0.01,
  label_cut_off = 0,
  label_col = NULL,
  ylim = c(0,1)
){
  # Format profile
  pf <- profile %>%
    separate(rname, c("pspecies", "poperon"), sep = "_", extra = "drop", remove = F) %>%
    mutate(
      umi_bin_size = as.integer(gsub(".*ubs=", "", qname)),
      qname = gsub(";.*", "", qname),
      poperon = gsub(";.*", "", poperon)
    ) %>%
    filter(grepl(species, rname))
  
  if (nrow(pf) > 0){
    # Count coverage of each individual operon
    pf_oc <- pf %>%
      distinct(qname, .keep_all = T) %>%
      group_by(pspecies, poperon) %>%
      summarise(
        operon_count = n(),
        ref_start  = min(ref_alnstart),
        ref_end = max(ref_alnlen + ref_alnstart - 1)
      )
    
    pf_oc_plot <- pf_oc %>%
      gather(key = pos,
             value = ref_pos,
             ref_start:ref_end) %>%
      mutate(err_freq = 0)
      
    
    # Calculate error frequency pr type and position
    pf_ec <- pf %>%
      group_by(pspecies, poperon, ref_pos, homopolymer, hpos, category, type) %>%
      summarise(err_count = n()) %>%
      left_join(pf_oc, by = c("pspecies", "poperon")) %>%
      mutate(
        err_freq = err_count/operon_count,
        hp_type = substr(homopolymer, 1,1)
      ) %>%
      filter(category != "perfect")
    
    
    # Baseplot
    p <- ggplot(pf_oc_plot,
                aes(
                  x = ref_pos,
                  y = err_freq,
                  )
                ) +
      # Add reference length
      geom_line(color = "darkgray") +
      geom_point(color = "darkgray") +
      # Add all data in gray
      geom_col(
        data = filter(
          pf_ec,
          err_count <= 1,
          err_freq <= lower_cut_off
          ),
        fill = "gray",
        color = "black",
        size = 0.1,
        width = 20,
        position = "identity"
        ) +
      # Add error data above defined threshold that is not singletons.
      geom_col(
        data = filter(
          pf_ec,
          err_count > 1,
          err_freq > lower_cut_off
          ),
        aes(fill = category),
        color = "black",
        size = 0.1,
        width = 20,
        position = "identity"
        ) +
      # Divide data accoring to taxa and operon
      facet_grid(rows = vars(poperon),
                 cols = vars(pspecies),
                 scales = "free",
                 space = "free_y",
                 drop = T) +
      # Set y-limits without removing data
      coord_cartesian(ylim = ylim) + 
      # Setup pretty plot labels
      theme_classic() +
      theme(text = element_text(size = 8),
            axis.text = element_text(size = 6)) +
      labs(fill = "Error type",
           x = "Reference position",
           y = "Error frequency")
    
    if(!is.null(label_col)){
      # Add annotation to data points i.e. type of error or specific position
      pf_ec_label <- filter(pf_ec, err_freq >= label_cut_off)
      p <- p +
        geom_text_repel(
          data = pf_ec_label,
          aes_string(
            label = label_col
          ),
          color = ifelse(is.na(pf_ec_label$hp_type), "black", "red"),
          direction = "x",
          segment.color = "grey50",
          segment.size = 0.2,
          nudge_y = ylim[2],
        )
    }
    
    # return plot
    return(p)
  } else {
    print(paste("No sequences of ",
                species,
                " with errors in current dataset.", sep = ""))
  }
}

# Process step error rates
lu_process_error_summary <- function(
  read_error,
  qc,
  ubs_cutoff,
  title,
  manual_error = rep(NA, 3),
  manual_sd = rep(NA, 3)
){
  # Filter data
  qc_f <- qc %>%
    filter(umi_bin_size >= ubs_cutoff)
  
  #Calculate statitics
  err_stats <- tibble(
    Step = c("Raw reads",
             paste("Consensus\n(bin size >=", ubs_cutoff,")", sep = ""),
             "Chimera rate"),
    `Mean\nerror (%)` = c(
      round(mean(read_error$ref_error),4),
      round(mean(qc_f$ref_error), 4),
      paste(sum(qc_f$chimera == "Y", na.rm = T),"/",length(qc_f$chimera), sep = "")
    ),
    `SD\nerror (%)` = c(
      round(sd(read_error$ref_error),4),
      round(sd(qc_f$ref_error), 4),
      ""
    )
  )
  
  # Add manual values if present
  err_stats <- err_stats %>%
    mutate( 
      `Mean\nerror (%)` = ifelse(
        is.na(manual_error),
        `Mean\nerror (%)`,
        manual_error
      ),
      `SD\nerror (%)` = ifelse(
        is.na(manual_sd),
        `SD\nerror (%)`,
        manual_sd
      )
    )
  
  # plot table
  grid.arrange(tableGrob(err_stats, rows= NULL, theme = tt), top = title)
}

# Plot errors as function of references
lu_ref_error_plot <- function(
  profile,
  title = NULL
){
  # Format data
  pf <- profile %>%
    ungroup() %>%
    mutate(
      bin_size = as.integer(gsub(".*=", "", qname)),
      qname = gsub(";.*", "", qname),
      species = gsub("_.*", "", rname),
      operon_n = gsub("^[a-zA-Z.0-9]*_|;.*|_.*", "", rname) %>%
        as.integer()
    ) %>%
    arrange(operon_n) %>%
    mutate(
      operon = as_factor(gsub("^[a-zA-Z.0-9]*_|;.*", "", rname))
    ) %>%
    spread(key = category, value = count) %>%
    mutate_at(vars(del:perfect), ~if_else(is.na(.), as.integer(0), .)) 
  
  #Summarise errors
  pfs <- pf %>%
    group_by(species, operon) %>% 
    summarise(nseqs = n(),
              mean_bin_size = round(mean(bin_size),0),
              median_bin_size = median(bin_size),
              median_errorc = median(ins+del+mm),
              mean_errorc = mean(ins+del+mm),
              median_mm = median(mm),
              mean_mm = mean(mm),
              median_ins = median(ins),
              mean_ins = mean(ins),
              median_del = median(del),
              mean_del = mean(del),       
              median_indel = median(ins+del),
              mean_indel = mean(ins+del)    
    )
  
  pf <- pf %>%
    group_by(qname, rname, species, operon) %>%
    summarise(n_err = sum(ins+del+mm))
  
  # Plot
  ggplot(
    pf,
    aes(
      x = operon,
      y = n_err,
      color = species,
      group = rname)
  ) +
    geom_jitter(
      width = 0.1,
      height = 0.1,
      alpha  = 0.05,
      size = 0.5
    ) +
    geom_errorbar(
      data = pfs, 
      aes(ymin = median_errorc, 
          ymax = median_errorc,
          x = operon), 
      color = "black",
      lwd = 0.25,
      width = 1,
      inherit.aes = F
    ) +
    facet_grid(~species, scales = "free_x", space = "free_x") +
    theme_classic() +
    labs(
      title = title,
      y = "Errors (n)",
      x = "Operon (#)"
      ) +
    theme(strip.background = element_blank(),
          legend.position = "none",
          text = element_text(size = 6),
          axis.text = element_text(size = 4),
          axis.text.x = element_text(angle=45),
          axis.ticks = element_line(size = 0.25, color = "darkgrey"),
          axis.ticks.length = unit(2, "pt"),
          axis.line = element_line(size = 0.25)) +
    annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=0.25)
}

lu_ref_abundance <- function(
  qc,
  title = NULL
){
  # Calculate relative abundance 
  abundance <- qc %>%
    group_by(ref_tax) %>%
    summarise(coverage = n()) %>%
    ungroup() %>%
    mutate(
      genus = gsub("_.*", "", ref_tax),
      operon = gsub("[A-Za-z]{1,}_|;.*|_.*", "", ref_tax) %>% as.integer(),
      operon_n = gsub(".*;size=|;", "", ref_tax) %>% as.integer(),
      coverage = coverage/operon_n
      ) %>%
    arrange(operon) %>%
    mutate(
      operon = as_factor(gsub("^[a-zA-Z.0-9]*_|;.*", "", ref_tax))
    )

  # Plot
  ggplot(
    abundance,
    aes(
      x = operon,
      y = coverage,
      color = genus
      )
    ) +
    geom_point(size = 0.5) +
    #scale_y_continuous(limits = c(0,110), breaks = c(0,25,50,75)) +
    facet_grid(~genus, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = 'off') +
    theme_classic() +
    labs(
      title = title,
      y = "Coverage (x)",
      x = "Operon (#)"
    ) +
    theme_classic() +
    theme(strip.background = element_blank(),
          legend.position = "none",
          text = element_text(size = 6),
          axis.text = element_text(size = 4),
          axis.text.x = element_text(angle=45),
          axis.ticks = element_line(size = 0.25, color = "darkgrey"),
          axis.ticks.length = unit(2, "pt"),
          axis.line = element_line(size = 0.25))
}

# Legacy code #################################################################
# 
# # Determine homopolymer positions in references
# lu_ref_hp_pos_old <- function(input){
#   refs_homo <- str_extract_all(string = input, pattern = c("A[A]+A|T[T]+T|C[C]+C|G[G]+G"))
#   names(refs_homo) <- paste(names(input), "x x")
#   position <- gregexpr(text = input, pattern = c("A[A]+A|T[T]+T|C[C]+C|G[G]+G"))
#   
#   homo <- data.frame(refname = names(unlist(refs_homo)), 
#                      homopolymer = unlist(refs_homo), 
#                      startpos = unlist(position),
#                      stringsAsFactors = F) %>%
#     separate(refname, into = "reference", sep = " x x", extra = "drop") %>%
#     mutate(hlen = paste("H",nchar(homopolymer), sep = ""),
#            hbase = substring(homopolymer, 1, 1),
#            stoppos = startpos + nchar(homopolymer)-1) %>%
#     select(reference, homopolymer, startpos, stoppos, hlen, hbase)
#   
#   return(homo)
# }
# 
# # Convert CS tag from sam file
# lu_parse_sam_cs_old <- function(input) {
#   d <- data.frame(qname = as.character(),                   ## Define output structure
#                   rname = as.character(), 
#                   flag = as.integer(), 
#                   alnlen = as.integer(),
#                   ref_alnlen = as.integer(),
#                   category = as.character(), 
#                   type = as.character(), 
#                   ref_pos = as.double(),
#                   ref_start_pos = as.integer())
#   
#   for (i in 1:nrow(input)){                                ## Loop through all sequences 
#     cs_in <- input$CS[i]
#     alnstart <- input$POS[i]
#     alnlen <- str_extract_all(input$CIGAR[i], "\\d+[MID]")[[1]] %>% 
#       str_extract( "(\\d+)") %>%
#       as.integer() %>%
#       sum()
#     ref_alnlen <- str_extract_all(input$CIGAR[i], "\\d+[MD]")[[1]] %>% 
#       str_extract( "(\\d+)") %>%
#       as.integer() %>%
#       sum()
#     
#     cs_in_hits <- gregexpr(text = cs_in, pattern = "\\-|\\+|\\*") %>% unlist()
#     
#     temp <- data.frame(cs = rep(cs_in, length(cs_in_hits)), start = cs_in_hits) %>%
#       mutate(type = substr(x = cs, start = start, stop = start),
#              category = ifelse(type == "*", "mm", ifelse(type == "-", "del", "ins")),
#              sub_cs = substr(x = cs, start = 0, stop = start-1),
#              type = ifelse(type == "*", substr(x = cs, start = start+1, stop = start+2), 
#                            ifelse(type == "-", substr(x = cs, start = start, stop = nchar(as.character(cs))) %>%
#                                     str_extract(pattern = c("\\-[a|t|c|g]+")), 
#                                   substr(x = cs, start = start, stop = nchar(as.character(cs))) %>%
#                                     str_extract(pattern = c("\\+[a|t|c|g]+"))))) %>%
#       rowwise() %>%
#       mutate(tdel = str_match_all(string = sub_cs, pattern = c("\\-[a|t|c|g]+")) %>% 
#                unlist() %>% 
#                str_count(pattern = "[a|t|c|g]") %>% 
#                sum(),
#              tmat = str_match_all(string = sub_cs, pattern = c("[0-9]+")) %>% 
#                unlist() %>% 
#                as.numeric() %>%
#                sum(),
#              tmm = str_count(sub_cs, pattern = "\\*"),
#              ref_pos = alnstart+tdel+tmat+tmm) %>%
#       mutate(qname = input$QNAME[i], 
#              rname = input$RNAME[i],
#              flag = input$FLAG[i],
#              alnlen = alnlen,
#              ref_alnstart = alnstart,
#              ref_alnlen = ref_alnlen) %>%
#       select(qname, rname, flag, alnlen, ref_alnstart, ref_alnlen, category, type, ref_pos) %>%
#       mutate(count = ifelse(category == "ins", nchar(type)-1, 1)) %>%
#       ungroup()
#     
#     d <- rbind.data.frame(d, temp)
#   }
#   
#   d$count[is.na(d$type)] <- 0
#   d$count <- as.integer(d$count)
#   d$category[is.na(d$type)] <- "perfect"
#   
#   #Inflate deltions to show at the positions they occur - this is done to characterise systematic errors better
#   temp2 <- d
#   
#   for (i in 1:nrow(d)){
#     if(d$category[i] == "del"){
#       if (nchar(d$type[i]) > 2){
#         for (j in 1:(nchar(d$type[i])-2)){
#           d[i,"ref_pos"] <-  d[i,"ref_pos"]+1
#           temp2 <- rbind.data.frame(temp2, d[i,])
#         }
#       }
#     } 
#   }
#   return(temp2)
# }
# 
# # Combine reference homopolymer positions with SNPs positions
# lu_error_profile_old <- function(refs, sam){
#   # Import data
#   sam_names = c("QNAME",
#                 "FLAG",
#                 "RNAME",
#                 "POS",
#                 "MAPQ",
#                 "CIGAR",
#                 "MNAME",
#                 "MPOS",
#                 "TLEN",
#                 "NM",
#                 "CS")
#   
#   if (is.character(sam)){
#     sam <- read_delim(file = sam,
#                       delim = "\t",
#                       col_names = sam_names)
#   }
#   
#   if (is.character(refs)){
#     refs <- readDNAStringSet(filepath = refs)
#   }
#   
#   # Detect homopolymer regions
#   homo <- lu_ref_hp_pos_old(refs)
#   
#   # Profile errors
#   sam_errors <- lu_parse_sam_cs_old(sam)
#   
#   d <- data.frame(sam_errors, homopolymer = "No", hpos = NA, stringsAsFactors = F)
#   
#   for (i in 1:nrow(sam_errors)){
#     temp <- filter(homo, reference == sam_errors$rname[i] & 
#                      startpos <= sam_errors$ref_pos[i] & 
#                      stoppos >= sam_errors$ref_pos[i])
#     if(nrow(temp)==1){d$homopolymer[i] <- temp$homopolymer}
#     if(nrow(temp)==1){d$hpos[i] <- temp$startpos}
#   }
#   return(d)
# } 
