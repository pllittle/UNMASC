# Code to visualize UNMASC data

rm(list=ls())
git_dir = "../../"
source(file.path(git_dir,"smart","RR","SOURCE.R"))
source(file.path(git_dir,"UNMASC","R","UNMASC.R"))

# ----------
# UNMASC definitions
# ----------
aa = UNMASC_definitions()
data.table::fwrite(aa,"../other/UNMASC_definitions.tsv",sep = "\t")

# ----------
# Simulate BAF and ARTI segments
# ----------
mix_AD = function(cc,DP,arti_maf,baf,psi){
	num_cnts = length(DP)
	
	if( cc == 1 ){
		AD = sapply(DP,function(xx) sample(seq(0,xx),1))
	} else if( cc == 2 ){
		if( psi == 0 ){
			AD = rbinom(n = num_cnts,size = DP,prob = 0.5)
		} else {
			AD = emdbook::rbetabinom(n = num_cnts,
				prob = 0.5,size = DP,theta = 1/psi)
		}
	} else if( cc == 3 ){
		if( psi == 0 ){
			AD = rbinom(n = num_cnts,size = DP,prob = baf)
		} else {
			AD = emdbook::rbetabinom(n = num_cnts,
				prob = baf,size = DP,theta = 1/psi)
		}
	} else if( cc == 4 ){
		if( psi == 0 ){
			AD = rbinom(n = num_cnts,size = DP,prob = 1 - baf)
		} else {
			AD = emdbook::rbetabinom(n = num_cnts,
				prob = 1 - baf,size = DP,theta = 1/psi)
		}
	} else if( cc == 5 ){
		if( psi == 0 ){
			AD = rbinom(n = num_cnts,size = DP,prob = arti_maf)
		} else {
			AD = emdbook::rbetabinom(n = num_cnts,
				prob = arti_maf,size = DP,theta = 1/psi)
		}
	}
	
	AD
}

num_segs = 8
mean_DP = 5e2
phi = 0.1
arti_maf = 0.05
psi = 0.01
bal_prob = c(2,1); bal_prob = bal_prob / sum(bal_prob)
h2m_prob = c(3,1); h2m_prob = h2m_prob / sum(h2m_prob)

set.seed(1); dat = c()
for(seg in seq(num_segs)){
	num_loci = sample(seq(5e2,2e3,500),1); num_loci
	seg_DP = rnbinom(n = num_loci,mu = mean_DP,size = 1/phi)
	
	seg_prob = c(0.1,0,0,0,0.1)
	balance = sample(c(FALSE,TRUE),1,prob = bal_prob)
	if( balance ){
		seg_prob[2] = 8
		seg_prob[3:4] = 0
		seg_baf = 0.5
	} else {
		seg_prob[2] = 0
		seg_prob[3:4] = 8
		seg_baf = runif(1,0.1,0.4)
	}
	baf_prob = seg_prob[2:4]
	baf_prob = baf_prob / sum(baf_prob)
	seg_prob[2:4] = baf_prob
	seg_prob = seg_prob / sum(seg_prob)
	seg_prob
	
	seg_cc = sample(x = seq(5),size = num_loci,
		replace = TRUE,prob = seg_prob)
	
	seg_AD = rep(NA,num_loci)
	seg_GERM = rep(0,num_loci)
	for(cc in seq(5)){
		idx = which(seg_cc == cc)
		seg_AD[idx] = mix_AD(cc = cc,DP = seg_DP[idx],
			arti_maf = arti_maf,baf = seg_baf,psi = psi)
		if( balance && cc == 2 ){
			seg_GERM[idx] = runif(length(idx),0.1,0.9)
		} else if( !balance && cc %in% c(3,4) ){
			seg_GERM[idx] = runif(length(idx),0.1,0.9)
		}
	}
	
	tmp_df = smart_df(seg = seg,cc = seg_cc,
		AD = seg_AD,DP = seg_DP)
	tmp_df$RD = seg_DP - seg_AD
	tmp_df$VAF = tmp_df$AD / tmp_df$DP
	tmp_df$GERM_AF = seg_GERM
	tmp_df$pt_cex = 0.5*log10(tmp_df$DP)
	dat = rbind(dat,tmp_df)
	rm(tmp_df,seg_cc,seg_GERM,seg_AD,seg_DP)
	
}

dim(dat); head(dat)

# hist(dat$DP,breaks = 40,col = "deepskyblue")

png("../other/sim_arti.png",units = 'px',
	height = 2000,width = 3000,res = 250,
	type = 'cairo',pointsize = 20)

mm = matrix(c(1,2,1,3),2,2,byrow = TRUE); mm
cex_lab = 1.2
layout(mm,width = c(2,1))
# layout.show(max(mm))
par(mar = c(4,4,1,0.5))
alpha = 0.3
dat$label_col = rgb(0,0,0,alpha)
dat$label_col[dat$GERM_AF > 0] = rgb(1,0,0,alpha)
dat$label_col[dat$cc == 5] = rgb(0,0,1,alpha)
plot(dat$VAF,cex = dat$pt_cex,pch = 16,
	bty = "n",col = dat$label_col,
	ylab = "Tumor VAF",cex.lab = cex_lab)
	abline(h = 0.5,lty = 2)
hist(dat$VAF,breaks = 50,col = "deepskyblue",
	xlab = "Tumor VAF",main = "All loci",cex.lab = cex_lab)
hist(dat$VAF[dat$GERM_AF == 0],breaks = 50,
	col = "deepskyblue",xlab = "Tumor VAF",
	main = "Non-germline DB loci",cex.lab = cex_lab)
par(mfrow = c(1,1),mar = c(5,4,4,2)+0.1)

dev.off()

###

