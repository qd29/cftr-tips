library(shiny)
library(bslib)
library(tools)
library(sangerseqR)
options(warn=-1)

# preview_text=t(HTML(paste("This is a preview version of <i>CFTR</i>-TIPS. Source codes are available <a href='https://github.com/qd29/cftr-tips' 
#                           target='_blank'>here</a> for full implementation."
#                           ,actionLink(inputId="preview_help",label=HTML("<span class='glyphicon glyphicon-question-sign'></span><br>")))))
preview_text=""

ui=fluidPage(
  theme=bs_theme(version=5),
  
  titlePanel(HTML("<i>CFTR</i>-TIPS: <i>CFTR</i> Tool for Inferring Poly-T/TG Size"),windowTitle="CFTR-TIPS"),
  preview_text,
  br(),
  h5(HTML(paste("Required information",actionLink(inputId="reqinfo_help",label=HTML("<span class='glyphicon glyphicon-question-sign'></span>"))))),
  fluidRow(
    column(width=4,fileInput(inputId="fwdab1","Forward ab1 file:")),
    column(width=4,fileInput(inputId="revab1","Reverse ab1 file:")),
  ),
  
  h5(HTML(paste("Optional information",actionLink(inputId="optinfo_help",label=HTML("<span class='glyphicon glyphicon-question-sign'></span>"))))),
  fluidRow(
    column(width=2,numericInput("mint","Minimum # of T:",value="3")),
    column(width=2,numericInput("maxt","Maximum # of T:",value="11")),
    column(width=2,numericInput("mintg","Minimum # of TG:",value="8")),
    column(width=2,numericInput("maxtg","Maximum # of TG:",value="16")),
  ),
  column(width=4,numericInput("minsig","Minimum informative Sanger trace signal:",value="200")),
  actionButton("run","Run Analysis",width="200px"),

  br(),
  t(HTML("<b>Demo data:</b>")), br(),
  t(HTML("Example 1 (5T-12TG/7T): <a href='https://github.com/qd29/cftr-tips/raw/main/demo_data/5T-12TG-7T_fwd.ab1'>forward trace</a>,
         <a href='https://github.com/qd29/cftr-tips/raw/main/demo_data/5T-12TG-7T_rev.ab1'> reverse trace</a>")), br(),
  t(HTML("Example 2 (5T-11TG/9T): <a href='https://github.com/qd29/cftr-tips/raw/main/demo_data/5T-11TG-9T_fwd.ab1'>forward trace</a>,
         <a href='https://github.com/qd29/cftr-tips/raw/main/demo_data/5T-11TG-9T_rev.ab1'> reverse trace</a>")), br(),
  t(HTML("Example 3 (5T-11TG/5T-12TG): <a href='https://github.com/qd29/cftr-tips/raw/main/demo_data/5T-11TG-5T-12TG_fwd.ab1'>forward trace</a>,
         <a href='https://github.com/qd29/cftr-tips/raw/main/demo_data/5T-11TG-5T-12TG_rev.ab1'> reverse trace</a>")), br(),

  plotOutput("img1",height="800px"),
  fluidRow(
    column(width=2,actionButton("prev","Previous",width="100%"),),
    column(width=2,actionButton("next1","Next",width="100%"),),
  ),
  br(),
  hr(),

  t(HTML("Contact: <a href='mailto:ding.qiliang@mayo.edu'>Qiliang (Andy) Ding, Ph.D.</a>")),
  br(),
  t(HTML("<b>© 2024 Mayo Foundation for Medical Education and Research</b> (Version: 20240316)")),
)

server=function(input,output){
  x=reactiveValues(rank=0)
  
  observeEvent(input$run,{
    x$rank=1
    
    fwdstart=NA
    fwdend=NA
    if(is.null(input$fwdab1$name)==F && toupper(file_ext(input$fwdab1$name))=="AB1"){
      fwdabi=sangerseq(read.abif(input$fwdab1$datapath))
      fwdsigtrace=fwdabi@traceMatrix
      fwdbasepos=fwdabi@peakPosMatrix[,1]
      
      fwdt=c()
      for (i in 1:(length(fwdbasepos)-1)){
        rt1=sum(fwdsigtrace[fwdbasepos[i]:fwdbasepos[i+1]-1,4])/sum(fwdsigtrace[fwdbasepos[i]:fwdbasepos[i+1]-1,])
        if (is.na(rt1)==F && rt1>=0.1){
          fwdt=c(fwdt,fwdbasepos[i])
        }
      }
      
      fwdseq=c(1,4,2,4,4,4,4,1,4,4,4,4,4,3,1) #ATCTTTTATTTTTGA
      for (i in 1:(length(fwdbasepos)-14)){
        fail=0
        for (j in 1:length(fwdseq)){
          rt1=fwdsigtrace[fwdbasepos[i+j-1],fwdseq[j]]/sum(fwdsigtrace[fwdbasepos[i+j-1],])
          if (is.na(rt1)==F && rt1<0.3){
            fail=1
          }
        }
        if (fail==0){
          fwdstart=i+15
        }
      }
      
      if (is.na(fwdstart)==F){
        for (i in fwdstart:length(fwdbasepos)){
          rt1=fwdsigtrace[fwdbasepos[i],4]/sum(fwdsigtrace[fwdbasepos[i],]) # modified from the Perl script
          if (rt1<0.3 && fwdbasepos[i]%in%fwdt==F){
            fwdend=i-1
            break
          }
        }
      }
      
      x$fwdsigtrace=fwdsigtrace
      x$fwdbasepos=fwdbasepos
      x$fwdstart=fwdstart
      x$fwdend=fwdend
    }
    
    revstart=NA
    revend=NA
    if(is.null(input$revab1$name)==F && toupper(file_ext(input$revab1$name))=="AB1"){
      revabi=sangerseq(read.abif(input$revab1$datapath))
      revsigtrace_orig=revabi@traceMatrix
      revsigtrace=cbind(revsigtrace_orig[,4],revsigtrace_orig[,3],revsigtrace_orig[,2],revsigtrace_orig[,1])
      revbasepos=revabi@peakPosMatrix[,1]
      
      revseq=c(3,3,3,3,4,4,4,1,3,3,3,1,2,1,1); #GGGGTTTAGGGACAA
      
      for (i in 1:(length(revbasepos)-14)){
        fail=0
        for (j in 1:length(revseq)){
          rt1=revsigtrace[revbasepos[i+j-1],revseq[j]]/sum(revsigtrace[revbasepos[i+j-1],])
          if (is.na(rt1)==F && rt1<0.3){
            fail=1
          }
        }
        if (fail==0){
          revstart=i+15
        }
      }
      
      if (is.na(revstart)==F){
        for (i in revstart:length(revbasepos)){
          rt1=revsigtrace[revbasepos[i],1]/sum(revsigtrace[revbasepos[i],])
          if (is.na(rt1)==F && rt1>=0.3){
            revend=i-1
            break
          }
        }
      }
      
      x$revsigtrace=revsigtrace
      x$revbasepos=revbasepos
      x$revstart=revstart
      x$revend=revend
    }

    if (is.na(fwdstart)==F && is.na(fwdend)==F && is.na(revstart)==F && is.na(revend)==F){
      indseqname=c()
      indseq=c()
      for (i in input$mintg:input$maxtg){
        for (j in input$mint:input$maxt){
          str1=paste0(j,"T-",i,"TG")
          str2=paste0(c(rep("TG",i),rep("T",j)),collapse="")
          indseqname=c(indseqname,str1)
          indseq=c(indseq,str2)
        }
      }
      
      combseqname=c()
      combseqfwd=c()
      combseqrev=c()
      for (i in 1:length(indseqname)){
        for (j in i:(length(indseqname)-1)){
          str1=""
          str2=paste(indseqname[i],indseqname[j],sep="/")
          arr1=strsplit(indseq[i],split="")[[1]]
          arr2=strsplit(indseq[j],split="")[[1]]
          maxlen=max(length(arr1),length(arr2))
          length(arr1)=maxlen
          length(arr2)=maxlen
          for (a in 1:maxlen){
            if (is.na(arr1[a])==F && is.na(arr2[a])==F){
              if (arr1[a]==arr2[a]){
                str1=paste0(str1,arr1[a])
              }
              else{
                str1=paste0(str1,"K")
              }
            }
            else{
              if (is.na(arr1[a])==F){
                str1=paste0(str1,arr1[a])
              }
              else{
                str1=paste0(str1,arr2[a])
              }
            }
          }
          
          str3=""
          arr1=strsplit(indseq[i],split="")[[1]]
          arr2=strsplit(indseq[j],split="")[[1]]
          arr3=rev(arr1)
          arr4=rev(arr2)
          maxlen=max(length(arr3),length(arr4))
          length(arr3)=maxlen
          length(arr4)=maxlen
          arr3=rev(arr3)
          arr4=rev(arr4)
          for (a in 1:maxlen){
            if (is.na(arr3[a])==F && is.na(arr4[a])==F){
              if (arr3[a]==arr4[a]){
                str3=paste0(str3,arr3[a])
              }
              else{
                str3=paste0(str3,"K")
              }
            }
            else{
              if (is.na(arr3[a])==F){
                str3=paste0(str3,arr3[a])
              }
              else{
                str3=paste0(str3,arr4[a])
              }
            }
          }
          combseqname=c(combseqname,str2)
          combseqfwd=c(combseqfwd,str1)
          combseqrev=c(combseqrev,str3)
        }
      }
      
      score=c()
      scorecomb=c()
      scoreseqfwd=c()
      scoreseqrev=c()
      for (i in 1:length(combseqname)){
        fail3=0
        ct=0
        tot=0
        arr5=strsplit(combseqfwd[i],split="")[[1]]
        length(arr5)=fwdend-fwdstart+1
        for (j in fwdstart:fwdend){
          if (fwdsigtrace[fwdbasepos[j],3]+fwdsigtrace[fwdbasepos[j],4]>=input$minsig){
            obs_g=fwdsigtrace[fwdbasepos[j],3]/(fwdsigtrace[fwdbasepos[j],3]+fwdsigtrace[fwdbasepos[j],4])
            obs_t=fwdsigtrace[fwdbasepos[j],4]/(fwdsigtrace[fwdbasepos[j],3]+fwdsigtrace[fwdbasepos[j],4])
            if (is.na(arr5[j-fwdstart+1])==F){
              if (arr5[j-fwdstart+1]=="T"){
                exp_g=0
                exp_t=1
                if (obs_t<0.3){
                  fail3=1
                }
              }
              else if (arr5[j-fwdstart+1]=="G"){
                exp_g=1
                exp_t=0
                if (obs_g<0.3){
                  fail3=1
                }
              }
              else{
                exp_g=0.5
                exp_t=0.5
                if (obs_t<0.1 || obs_g<0.1){
                  fail3=1
                }
              }
              ct=ct+2
              tot=tot+(obs_g-exp_g)^2+(obs_t-exp_t)^2
            }
          }
        }
        
        arr6=strsplit(combseqrev[i],split="")[[1]]
        arr6=rev(arr6)
        length(arr6)=revend-revstart+1
        for (j in revstart:revend){
          if (revsigtrace[revbasepos[j],3]+revsigtrace[revbasepos[j],4]>=input$minsig){
            obs_g=revsigtrace[revbasepos[j],3]/(revsigtrace[revbasepos[j],3]+revsigtrace[revbasepos[j],4])
            obs_t=revsigtrace[revbasepos[j],4]/(revsigtrace[revbasepos[j],3]+revsigtrace[revbasepos[j],4])
            if (is.na(arr6[j-revstart+1])==F){
              if (arr6[j-revstart+1]=="T"){
                exp_g=0
                exp_t=1
                if (obs_t<0.3){
                  fail3=1
                }
              }
              else if (arr6[j-revstart+1]=="G"){
                exp_g=1
                exp_t=0
                if (obs_g<0.3){
                  fail3=1
                }
              }
              else{
                exp_g=0.5
                exp_t=0.5
                if (obs_t<0.1 || obs_g<0.1){
                  fail3=1
                }
              }
              ct=ct+2
              tot=tot+(obs_g-exp_g)^2+(obs_t-exp_t)^2
            }
          }
        }
        
        if (fail3==0 && nchar(combseqfwd[i])<=fwdend-fwdstart+1 && nchar(combseqrev[i])>=revend-revstart+1){
          score=c(score,round(tot*1e3/ct)/1e3+(100-nchar(combseqfwd[i])-nchar(combseqrev[i]))/100000)
          scorecomb=c(scorecomb,combseqname[i])
          scoreseqfwd=c(scoreseqfwd,combseqfwd[i])
          scoreseqrev=c(scoreseqrev,combseqrev[i])
        }
      }
      
      if(length(score)>=1){
        scoreindex=sort(score,index=T)$ix
        x$score=score
        x$scorecomb=scorecomb
        x$scoreindex=sort(score,index=T)$ix
        x$scoreseqfwd=scoreseqfwd
        x$scoreseqrev=scoreseqrev
        x$rank=1
      }
      else{
        x$rank=10001
      }
    }
    else{
      x$rank=10002
    }
   
  })
    
  observeEvent(input$prev,{
    if(x$rank>1 && x$rank<=10000){
      x$rank=x$rank-1
    }
  })
  
  observeEvent(input$next1,{
    if(x$rank>=1 && x$rank<=10000 && x$rank<length(x$scoreindex)){
      x$rank=x$rank+1
    }
  })
  
  output$img1=renderPlot({
    if(x$rank==0){
      par(mfrow=c(3,1),mar=c(2,4,2,2))
      plot(0,0,col=rgb(0,0,0,0),axes=F,xlab="",ylab="",xlim=c(0,10),ylim=c(0,10),xaxs="i")
      text(5,5,"Please upload the forward and reverse ab1 files,\n press the Run Analysis button,\n and wait for the results to load.",cex=4)
      plot(0,0,col=rgb(0,0,0,0),axes=F,xlab="",ylab="",xlim=c(0,10),ylim=c(0,10),xaxs="i")
      plot(0,0,col=rgb(0,0,0,0),axes=F,xlab="",ylab="",xlim=c(0,10),ylim=c(0,10),xaxs="i")
    }
    else if(x$rank==10001){
      par(mfrow=c(3,1),mar=c(2,4,2,2))
      plot(0,0,col=rgb(0,0,0,0),axes=F,xlab="",ylab="",xlim=c(0,10),ylim=c(0,10),xaxs="i")
      text(5,5,"No combinations match the uploaded Sanger traces.\nConsider increase the range of T and TG repeats.",cex=4)
      plot(0,0,col=rgb(0,0,0,0),axes=F,xlab="",ylab="",xlim=c(0,10),ylim=c(0,10),xaxs="i")
      plot(0,0,col=rgb(0,0,0,0),axes=F,xlab="",ylab="",xlim=c(0,10),ylim=c(0,10),xaxs="i")
    }
    else if(x$rank==10002){
      par(mfrow=c(3,1),mar=c(2,4,2,2))
      plot(0,0,col=rgb(0,0,0,0),axes=F,xlab="",ylab="",xlim=c(0,10),ylim=c(0,10),xaxs="i")
      text(5,5,"Unable to detect the poly-T/TG tract\nin the forward and/or reverse trace.\n\nCheck if correct files were uploaded.\nFile names must end in .ab1.\nIf issue presists, consider manual review of the Sanger traces.",cex=4)
      plot(0,0,col=rgb(0,0,0,0),axes=F,xlab="",ylab="",xlim=c(0,10),ylim=c(0,10),xaxs="i")
      plot(0,0,col=rgb(0,0,0,0),axes=F,xlab="",ylab="",xlim=c(0,10),ylim=c(0,10),xaxs="i")
    }
    else if(x$rank>=1){
      par(mfrow=c(3,1),mar=c(2,4,2,2))
      arr7=strsplit(x$scorecomb[x$scoreindex[x$rank]],split="/")[[1]]
      arr7=sort(arr7)
      plot(0,0,col=rgb(0,0,0,0),axes=F,xlab="",ylab="",xlim=c(0,10),ylim=c(0,10),xaxs="i")
      text(0,9,paste0("Rank: ",x$rank,"     ",arr7[1],"/",arr7[2],"     Normalized Difference: ",
                      x$score[x$scoreindex[x$rank]]),cex=3.5,pos=4)
      text(0,6.5,paste0("Forward Trace: ",isolate(input$fwdab1$name)),cex=2,pos=4)
      text(0,5.2,paste0("Reverse Trace: ",isolate(input$revab1$name)),cex=2,pos=4)
      text(0,3.9,date(),cex=2,pos=4)
      
      fwdpolypos=c()
      fwdmaxlen=0
      for (i in 1:length(x$scoreseqfwd)){
        if (fwdmaxlen<nchar(x$scoreseqfwd[i])){
          fwdmaxlen=nchar(x$scoreseqfwd[i])
        }
      }
      
      for (i in 1:fwdmaxlen){
        arr8=c()
        for (j in 1:length(x$scoreseqfwd)){
          arr9=strsplit(x$scoreseqfwd[j],split="")[[1]]
          length(arr9)=fwdmaxlen
          arr8=c(arr8,arr9[i])
        }
        fwdpolypos=c(fwdpolypos,length(unique(arr8)))
      }
      
      plot(0,0,col=rgb(0,0,0,0),axes=F,xlab="",ylab="",xaxs="i",xlim=c(x$fwdbasepos[x$fwdstart-6],x$fwdbasepos[x$fwdend+6]),
           ylim=c(-max(x$fwdsigtrace[x$fwdbasepos[x$fwdstart-6]:x$fwdbasepos[x$fwdend+6]])*0.2,
                  max(x$fwdsigtrace[x$fwdbasepos[x$fwdstart-6]:x$fwdbasepos[x$fwdend+6]])*1.3),main="Forward Trace",cex.main=2)
      rect(x$fwdbasepos[x$fwdstart],-1000,x$fwdbasepos[x$fwdend],max(x$fwdsigtrace)*2,col="#f0f0f0",border=F)
      box()
      axis(2,cex.axis=1.5)
      lines(x$fwdsigtrace[,1],col="green")
      lines(x$fwdsigtrace[,2],col="blue")
      lines(x$fwdsigtrace[,3],lwd=2)
      lines(x$fwdsigtrace[,4],col="red",lwd=2)
      strfs=strsplit(x$scoreseqfwd[x$scoreindex[x$rank]],split="")[[1]]
      strfcol=rep("gray",length(strfs))
      strfcol[which(strfs=="T")]="red"
      strfcol[which(strfs=="G")]="black"
      strfcol[which(strfs=="K")]="darkred"
      strfs[which(strfs=="K")]="T/G"
      for(i in 1:fwdmaxlen){
        if (fwdpolypos[i]>1){
          rect(x$fwdbasepos[x$fwdstart+i-1]-5,-max(x$fwdsigtrace[x$fwdbasepos[x$fwdstart-6]:x$fwdbasepos[x$fwdend+6]])*0.2,
               x$fwdbasepos[x$fwdstart+i-1]+5,-max(x$fwdsigtrace[x$fwdbasepos[x$fwdstart-6]:x$fwdbasepos[x$fwdend+6]])*0.1,
               border=NA,col="#9ecae1")
        }
      }
      for (i in 1:length(strfs)){
        text(x$fwdbasepos[x$fwdstart+i-1],-max(x$fwdsigtrace[x$fwdbasepos[x$fwdstart-6]:x$fwdbasepos[x$fwdend+6]])*0.15,label=strfs[i],col=strfcol[i],cex=2)
      }
     
      revpolypos=c()
      revmaxlen=0
      for (i in 1:length(x$scoreseqrev)){
        if (revmaxlen<nchar(x$scoreseqrev[i])){
          revmaxlen=nchar(x$scoreseqrev[i])
        }
      }
      
      for (i in 1:revmaxlen){
        arr8=c()
        for (j in 1:length(x$scoreseqrev)){
          arr9=rev(strsplit(x$scoreseqrev[j],split="")[[1]])
          length(arr9)=revmaxlen
          arr8=c(arr8,arr9[i])
        }
        revpolypos=c(revpolypos,length(unique(arr8)))
      }
      
      plot(0,0,col=rgb(0,0,0,0),axes=F,xlab='',ylab='',xaxs='i',xlim=c(x$revbasepos[x$revend+6],x$revbasepos[x$revstart-6]),
           ylim=c(-max(x$revsigtrace[x$revbasepos[x$revend+6]:x$revbasepos[x$revstart-6]])*0.2,
                  max(x$revsigtrace[x$revbasepos[x$revend+6]:x$revbasepos[x$revstart-6]])*1.3),main='Reverse Trace',cex.main=2)
      rect(x$revbasepos[x$revstart],-1000,x$revbasepos[x$revend],max(x$revsigtrace)*2,col="#f0f0f0",border=F)
      box()
      axis(2,cex.axis=1.5)
      lines(x$revsigtrace[,1],col="green")
      lines(x$revsigtrace[,2],col="blue")
      lines(x$revsigtrace[,3],lwd=2)
      lines(x$revsigtrace[,4],col="red",lwd=2)
      strrs=strsplit(x$scoreseqrev[x$scoreindex[x$rank]],split="")[[1]]
      strrs=rev(strrs)
      strrcol=rep("gray",length(strrs))
      strrcol[which(strrs=="T")]="red"
      strrcol[which(strrs=="G")]="black"
      strrcol[which(strrs=="K")]="darkred"
      strrs[which(strrs=="K")]="T/G"
      for(i in 1:revmaxlen){
        if (revpolypos[i]>1){
          rect(x$revbasepos[x$revstart+i-1]-5,-max(x$revsigtrace[x$revbasepos[x$revend+6]:x$revbasepos[x$revstart-6]])*0.2,
               x$revbasepos[x$revstart+i-1]+5,-max(x$revsigtrace[x$revbasepos[x$revend+6]:x$revbasepos[x$revstart-6]])*0.1,
               border=NA,col="#9ecae1")
        }
      }
      for (i in 1:length(strrs)){
        text(x$revbasepos[x$revstart+i-1],-max(x$revsigtrace[x$revbasepos[x$revend+6]:x$revbasepos[x$revstart-6]])*0.15,label=strrs[i],col=strrcol[i],cex=2)
      }
    }
  })
  
  observeEvent(input$reqinfo_help,{
    showModal(modalDialog(title="Help - Required Information",size="l",easyClose=T,
                          HTML("Please upload the forward and reverse Sanger chromatograms of the same sample.<br>
                               At a minimum, the chromatograms must cover the (TG)mTn tract ± 15 bp.<br>
                               The file extensions <b>must</b> be in .ab1.")))
  })
  
  observeEvent(input$optinfo_help,{
    showModal(modalDialog(title="Help - Optional Information",size="l",easyClose=T,
                          HTML("The search space of possible (TG)mTn allele combinations can be re-defined using the <b>“Minimum # of T”, “Maximum # of T”,
                          “Minimum # of TG”, and “Maximum # of TG”</b> parameters. Given that the default search space encompasses all known (TG)mTn alleles, 
                               adjustments to these parameters will be rarely necessary, if at all.<br><br>
                               You may also adjust the <b>“Minimum informative Sanger trace signal”</b> parameter. In the Sanger chromatograms, positions 
                               with signal intensity below this value will not be used to compare the observed vs. expected patterns. We recommend adjusting 
                               this parameter based on the overall quality of your Sanger chromatograms.")))
  })
  
  observeEvent(input$preview_help,{
    showModal(modalDialog(title="Limitations - Preview Version",size="l",easyClose=T,
                          HTML("1. <i>CFTR</i>-TIPS responds slower in the preview version than the full implementation on local computers. Occasionally, the
                          preview website may become unresponsive. In this scenario, we recommend trying again in a few minutes. These issues are due to
                          resource limitations of the website hosting service.<br><br>
                          2. In the preview version, you must reload the webpage if you would like to analyze chromatograms of another sample. In other words,
                          the preview website does not allow repeated upload of chromatograms. The full implementation does not have this limitation.
                          <br><br>3. While we do not retain any uploaded data, de-identifying your samples when using the preview version is recommended.")))
  })
}

runApp(list(ui=ui,server=server),launch.browser = TRUE) # for local implementations
# shinyApp(ui,server) # for shinyapps.io