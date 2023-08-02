          seed = -1
       seqfile = pres90.aa.noog.phy
      treefile = aa_for_mcmc_time.tre
      mcmcfile = culicidmcmc.txt
       outfile = culicid.txt

         ndata = 1
       seqtype = 2    * 0: nucleotides; 1:codons; 2:AAs
       usedata = 2    * 0: no data; 1:seq like; 2:normal approximation; 3:out.BV (in.BV)
         clock = 2    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge = '<250'  * safe constraint on root age, used if no fossil for root.

         model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0    * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

         cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

         BDparas = 1 1 0   * birth1, death1, sampling0
         kappa_gamma = 6 2     * gamma prior for kappa
         alpha_gamma = 1 1     * gamma prior for alpha

         rgene_gamma = 2 40 1   * gammaDir prior for rate for genes
         sigma2_gamma = 1 10 1   * gammaDir prior for sigma^2     (for clock=2 or 3)

         print = 1   * 0: no mcmc sample; 1: everything except branch rates 2: everything
         burnin = 200000
         sampfreq = 4000
         nsample = 10000
