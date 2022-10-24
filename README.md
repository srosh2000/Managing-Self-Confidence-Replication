# Managing-Self-Confidence-Replication


This describes the replication package for "The Role of Gender Heterogeneity and Social Learning in Belief Updating" by Roshini Sudhaharan. This is partly a replication study of Mobius et al.,(2022).

All figures, tables and statistics reported in the manuscript are
produced from two data files:

    /data/main_session.csv          Contains data from the main and followup experiments on belief updating
    /data/competition_session.csv   Contains data from the second follow-up experiment on competitive behavior

The variables in these datasets are described in corresponding
codebook files (/data/*_codebook.xlsx), and a copy of the
experimental instructions is also provided
(/data/experimental_instructions.pdf). Information about subject
eligibility, recruitment and demographics is available in the main
paper itself (Section 3).

Results are produced from these datasets using two scripts:

    /code/analysis_main_RS.r           Produces the bulk of the results in the paper.
    /code/analysis_competition_RS.r    Produces some additional results that appear in the section  on competitive behavior.
    

Within each code file, sections are demarcated with comments that
indicate what piece of paper content they produce (e.g. "Figure
1"). The scripts make use of publicly available R libraries and
also of three user-defined libraries (lib_clmclx.r, lib_hacktex.r,
and lib_multiregtable.r) which are also included in /code/. Figures
and tables produced by the scripts are stored in /figures and
/tables, respectively.
