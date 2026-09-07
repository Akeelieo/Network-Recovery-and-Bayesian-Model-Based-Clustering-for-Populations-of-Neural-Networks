#-------------------------------------------------------------------------------
# EEG preparation functions for the stochastic multi-population JRNMM
#
# Replaces the hard-coded prepare_EEG_data1.R / prepare_EEG_data2.R scripts.
# Seizure times are read from the CHB-MIT summary file (e.g. chb01-summary.txt)
# instead of being typed in by hand.
#
# The numerical preprocessing is IDENTICAL to prepare_EEG_data1.R:
#   - the same sample indices          l = floor(t_start/hreal), r = floor(t_end/hreal)-1
#   - the same scaling                 sc = 0.05
#   - the same interpolation           approx() onto gridsim, then na_interpolation()
# so results remain comparable with runs produced by the original scripts.
#-------------------------------------------------------------------------------

library(edf)
library(imputeTS)

#-------------------------------------------------------------------------------
# 1. Parse a CHB-MIT summary file
#-------------------------------------------------------------------------------
#
# Returns a data.frame with one row per seizure:
#   file            e.g. "chb01_03.edf"
#   record          e.g. "chb01_03"
#   seizure_index   1, 2, ... within that file
#   start, end      seizure onset/offset in seconds from the start of the file
#   duration        end - start
#
# Handles both label styles found in the database:
#   "Seizure Start Time: 2996 seconds"
#   "Seizure 1 Start Time: 2996 seconds"

read_seizure_summary<-function(summary_path){

  if(!file.exists(summary_path)){
    stop("Summary file not found: ",summary_path)
  }

  lines<-readLines(summary_path,warn=FALSE)
  lines<-trimws(lines)

  current_file<-NA_character_
  starts<-c(); ends<-c(); files<-c()

  for(ln in lines){

    if(grepl("^File Name:",ln)){
      current_file<-trimws(sub("^File Name:","",ln))
      next
    }

    #"Seizure Start Time:" or "Seizure 1 Start Time:"
    if(grepl("^Seizure( [0-9]+)? Start Time:",ln)){
      v<-as.numeric(gsub("[^0-9.]","",sub(".*Start Time:","",ln)))
      starts<-c(starts,v); files<-c(files,current_file)
      next
    }

    if(grepl("^Seizure( [0-9]+)? End Time:",ln)){
      v<-as.numeric(gsub("[^0-9.]","",sub(".*End Time:","",ln)))
      ends<-c(ends,v)
      next
    }
  }

  if(length(starts)==0){
    stop("No seizures found in ",summary_path)
  }
  if(length(starts)!=length(ends)){
    stop("Mismatched seizure start/end entries in ",summary_path,
         " (",length(starts)," starts, ",length(ends)," ends)")
  }

  out<-data.frame(file=files,
                  record=sub("\\.edf$","",files),
                  start=starts,
                  end=ends,
                  duration=ends-starts,
                  stringsAsFactors=FALSE)

  #index seizures within each file
  out$seizure_index<-ave(seq_len(nrow(out)),out$file,FUN=seq_along)
  out<-out[,c("file","record","seizure_index","start","end","duration")]
  rownames(out)<-NULL
  return(out)
}

#-------------------------------------------------------------------------------
# 2. Look up the window for one record / period
#-------------------------------------------------------------------------------
#
# period = "during" -> [onset, onset + T_window)
# period = "before" -> [onset - T_window, onset)
#
# Errors if the requested window does not fit: a seizure shorter than T_window,
# a "before" window running off the start of the file, or a "before" window
# overlapping the previous seizure in the same file.

get_window<-function(seizures,record,period,T_window=40,seizure_index=1){

  period<-match.arg(period,c("during","before"))

  row<-seizures[seizures$record==record & seizures$seizure_index==seizure_index,]
  if(nrow(row)==0){
    stop("No seizure ",seizure_index," found for record ",record,
         ". Available: ",paste(unique(seizures$record),collapse=", "))
  }

  onset<-row$start[1]; offset<-row$end[1]; dur<-row$duration[1]

  if(period=="during"){
    if(dur<T_window){
      warning("Seizure in ",record," lasts ",dur,"s, shorter than T_window=",T_window,
              "s. The window extends ",T_window-dur,"s past seizure end.",
              call.=FALSE, immediate.=TRUE)
    }
    if(dur>T_window){
      message("Note: seizure in ",record," lasts ",dur,
              "s; using only the first ",T_window,"s from onset.")
    }
    win_start<-onset
    } else {
    win_start<-onset-T_window
    if(win_start<0){
      stop("Before-seizure window for ",record," starts at ",win_start,
           "s, before the start of the file.")
    }
    #check it does not overlap an earlier seizure in the same file
    prev<-seizures[seizures$record==record & seizures$seizure_index<seizure_index,]
    if(nrow(prev)>0 && max(prev$end)>win_start){
      stop("Before-seizure window for ",record," overlaps the preceding seizure ",
           "(ends at ",max(prev$end),"s).")
    }
  }

  list(record=record,period=period,seizure_index=seizure_index,
       win_start=win_start,win_end=win_start+T_window,
       onset=onset,offset=offset,duration=dur,T_window=T_window)
}

#-------------------------------------------------------------------------------
# 3. Extract, scale and interpolate one window from an .edf file
#-------------------------------------------------------------------------------
#
# Returns an N x length(gridsim) matrix.
# channels are given with hyphens ("FP1-F7"); the edf package stores them with
# underscores, so they are converted internally.

extract_EEG_window<-function(edf_path,win,
                             channels=c("FP1-F7","FP1-F3","FP2-F4","FP2-F8"),
                             hsim=2*10^-3,hz=256,sc=0.05,edf_data=NULL){

  if(is.null(edf_data)){
    if(!file.exists(edf_path)) stop("EDF file not found: ",edf_path)
    message("Reading ",edf_path," (this can take a while) ...")
    edf_data<-read.edf(edf_path)
  }

  T_window<-win$T_window
  hreal<-1/hz

  gridsim <-seq(from=0,to=T_window,by=hsim)
  gridreal<-seq(from=0,to=T_window-hreal,by=hreal)

  #sample indices - identical convention to prepare_EEG_data1.R
  l<-floor(win$win_start/hreal)
  r<-floor(win$win_end/hreal)-1

  ch_names<-gsub("-","_",channels)
  available<-names(edf_data$signal)

  missing<-setdiff(ch_names,available)
  if(length(missing)>0){
    stop("Channel(s) not found in ",basename(edf_path),": ",
         paste(gsub("_","-",missing),collapse=", "),
         "\nAvailable: ",paste(gsub("_","-",available),collapse=", "))
  }

  N<-length(ch_names)
  X<-matrix(0,nrow=N,ncol=length(gridsim))

  for(i in 1:N){
    sig<-edf_data$signal[[ch_names[i]]]$data

    if(r>length(sig)){
      stop("Requested window runs past the end of ",basename(edf_path),
           " (need sample ",r,", file has ",length(sig),").")
    }

    seg<-sig[l:r]*sc

    if(length(seg)!=length(gridreal)){
      stop("Extracted segment for ",channels[i]," has length ",length(seg),
           " but the 256 Hz grid has length ",length(gridreal),".")
    }

    xi<-approx(x=gridreal,y=seg,xout=gridsim)$y   #onto the simulation grid
    xi<-na_interpolation(xi,option="linear")      #fill the trailing NA
    X[i,]<-xi
  }

  #--- validation
  if(any(!is.finite(X))) stop("Non-finite values in the extracted data.")
  sds<-apply(X,1,sd)
  if(any(sds==0)) stop("Channel(s) with zero variance: ",
                       paste(channels[sds==0],collapse=", "))
  message("Extracted ",N," channels x ",ncol(X),
          " points; sd = ",paste(round(sds,3),collapse=", "))

  attr(X,"edf_data")<-NULL
  return(X)
}

#-------------------------------------------------------------------------------
# 4. Wrapper: prepare (or load cached) reference data for one record/period
#-------------------------------------------------------------------------------
#
# Writes  <out_dir>/X1.txt ... XN.txt  and  <out_dir>/metadata.txt
# and returns a list with the data matrix and the paths.
#
# If the files already exist and force=FALSE, they are loaded instead of
# being regenerated (reading an .edf takes far longer than the inference does).

prepare_EEG_data<-function(record,period,
                           summary_path,edf_dir=".",out_root="RefData",
                           T_window=40,hsim=2*10^-3,hz=256,sc=0.05,
                           channels=c("FP1-F7","FP1-F3","FP2-F4","FP2-F8"),
                           seizure_index=1,force=FALSE){

  period<-match.arg(period,c("during","before"))
  period_dir<-ifelse(period=="during","duringSeizure","beforeSeizure")

  out_dir<-file.path(out_root,record,period_dir,"Data")
  gridsim<-seq(from=0,to=T_window,by=hsim)
  N<-length(channels)

  xfiles<-file.path(out_dir,paste0("X",1:N,".txt"))

  #--- use cached files if present
  if(!force && all(file.exists(xfiles))){
    message("Using cached data in ",out_dir)
    X<-matrix(0,nrow=N,ncol=length(gridsim))
    for(i in 1:N){
      X[i,]<-as.vector(t(read.table(xfiles[i],header=FALSE)))
    }
    return(list(X=X,dir=out_dir,gridsim=gridsim,cached=TRUE))
  }

  #--- otherwise generate
  seizures<-read_seizure_summary(summary_path)
  win<-get_window(seizures,record,period,T_window=T_window,
                  seizure_index=seizure_index)

  edf_path<-file.path(edf_dir,paste0(record,".edf"))
  X<-extract_EEG_window(edf_path,win,channels=channels,
                        hsim=hsim,hz=hz,sc=sc)

  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  LX<-ncol(X)
  for(i in 1:N){
    write(t(X[i,]),file=xfiles[i],ncolumns=LX,sep=" ")
  }

  #--- provenance
  meta<-c(
    paste0("record          : ",record),
    paste0("period          : ",period),
    paste0("seizure_index   : ",win$seizure_index),
    paste0("seizure_onset_s : ",win$onset),
    paste0("seizure_offset_s: ",win$offset),
    paste0("seizure_dur_s   : ",win$duration),
    paste0("window_start_s  : ",win$win_start),
    paste0("window_end_s    : ",win$win_end),
    paste0("T_window_s      : ",T_window),
    paste0("hsim            : ",hsim),
    paste0("sampling_hz     : ",hz),
    paste0("scale_factor    : ",sc),
    paste0("channels        : ",paste(channels,collapse=", ")),
    paste0("n_points        : ",LX),
    paste0("channel_sd      : ",paste(round(apply(X,1,sd),4),collapse=", ")),
    paste0("source_edf      : ",edf_path),
    paste0("summary_file    : ",summary_path),
    paste0("generated       : ",format(Sys.time(),"%Y-%m-%d %H:%M:%S"))
  )
  writeLines(meta,file.path(out_dir,"metadata.txt"))

  message("Wrote ",N," channel files and metadata.txt to ",out_dir)
  return(list(X=X,dir=out_dir,gridsim=gridsim,cached=FALSE,window=win))
}

#-------------------------------------------------------------------------------
# 5. Optional: plot the extracted reference data (as in prepare_EEG_data1.R)
#-------------------------------------------------------------------------------

plot_reference_data<-function(X,gridsim,out_dir,
                              channels=c("FP1-F7","FP1-F3","FP2-F4","FP2-F8"),
                              ylim_val=c(-20,20)){
  N<-nrow(X)
  pdf(width=8,height=4,file.path(out_dir,"referenceData.pdf"))
  par(mfrow=c(N+1,1),mai=c(0.01,0.65,0.05,0.3))
  for(j in 1:N){
    plot(gridsim,X[j,],type="l",col="blue",xlab="",ylab=channels[j],
         ylim=ylim_val,yaxt="n",xaxt="n")
    axis(2,c(-15,0,15))
  }
  axis(1,pretty(gridsim))
  dev.off()
  message("Wrote ",file.path(out_dir,"referenceData.pdf"))
}
