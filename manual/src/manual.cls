\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{manual}[2020/10/10]

%% Load base class
\LoadClass{report}

%% Load packages
\usepackage[toc,page]{appendix}
\usepackage{geometry}
\usepackage{import}
\usepackage{titlesec}
\usepackage{cmap}
\usepackage{amsmath,amsfonts,amssymb}
\usepackage{varwidth}
\usepackage{graphicx}
\usepackage{setspace}
\usepackage[version=3]{mhchem}
\usepackage{graphicx}
\usepackage{latexsym}
\usepackage{keyval}
\usepackage{ifthen}
\usepackage{moreverb}
\usepackage{xcolor}
\usepackage{multicol}
\usepackage{svg}
\usepackage{pdfpages} 
\usepackage{physics} 
\usepackage{caption}
\usepackage{subcaption}
\usepackage{lscape}
\usepackage{placeins}


% Equation packages
\usepackage{amsmath,amsfonts,amssymb}
\usepackage{textgreek}
\usepackage[nice]{nicefrac}

% Frames
\usepackage{mdframed}

%% Bib and Referencing
\usepackage[backend=bibtex,
bibstyle=numeric,
hyperref=true,
sorting=none]{biblatex}

\bibliography{manual.bib}  

\usepackage{notoccite}



%% Page settings
\pagenumbering{arabic}
\geometry{a4paper,total={210mm,297mm},left=20mm,right=20mm,top=20mm,bottom=20mm}

%% Margin settings
\def\changemargin#1#2{\list{}{\rightmargin#2\leftmargin#1}\item[]}
\let\endchangemargin=\endlist

%% Paragraph settings
\usepackage{setspace}
\setlength{\parindent}{0pt}
\setlength{\parskip}{1.3ex plus 0.5ex minus 0.3ex}
\onehalfspacing
%%\doublespacing

%% Tikz package
\usepackage{tikz}
\usetikzlibrary{shapes,arrows}

%% Colours
\usepackage{color}
\definecolor{red}{rgb}{1.0,0,0}
\definecolor{green}{rgb}{0,1.0,0}
\definecolor{blue}{rgb}{0,0,1.0}
\definecolor{dkgreen}{rgb}{0,0.6,0}
\definecolor{gray}{rgb}{0.5,0.5,0.5}
\definecolor{mauve}{rgb}{0.58,0,0.82}

%% Listings for code
\usepackage{listings}
\lstloadlanguages{Fortran,Bash}
\lstset{frame=tb,
  aboveskip=3mm,
  belowskip=3mm,
  showstringspaces=false,
  columns=flexible,
  basicstyle={\small\ttfamily}
}
%% Listing styles
\lstdefinestyle{sBash}
{
  numbers=left,
  stepnumber=1,
  language=bash,
  tabsize=4,
  keywordstyle=\color{red}
}
\lstdefinestyle{sFortran}
{
  numbers=left,
  stepnumber=1,
  language=Fortran,
  tabsize=4,
  numberstyle=\tiny\color{gray},
  keywordstyle=\color{blue},
  commentstyle=\color{dkgreen},
  stringstyle=\color{mauve}
}

%% Algorithm + Pseduocode
\usepackage{algorithmicx}
\usepackage{algpseudocode}
\algdef{SE}[DOWHILE]{Do}{doWhile}{\algorithmicdo}[1]{\algorithmicwhile\ #1}%


%%------------------------------------------------------------------------------
%%------------------------------------------------------------------------------
%% Define new commands

%% Thesis
%%-----------------------------
\def\changemargin#1#2{\list{}{\rightmargin#2\leftmargin#1}\item[]}
\let\endchangemargin=\endlist
\newcommand{\abstractpreamble}[1]{\textit{#1}}
\newcommand{\subheadings}{\noindent\textbf}





\usepackage{color}
\definecolor{red}{rgb}{1.0,0,0}
\definecolor{green}{rgb}{0,1.0,0}
\definecolor{blue}{rgb}{0,0,1.0}
\definecolor{dkgreen}{rgb}{0,0.6,0}
\definecolor{gray}{rgb}{0.5,0.5,0.5}
\definecolor{grey}{rgb}{0.5,0.5,0.5}
\definecolor{light-gray}{rgb}{0.9,0.9,0.9}
\definecolor{mauve}{rgb}{0.58,0,0.82}
\definecolor{white}{rgb}{1.0,1.0,1.0}


\definecolor{blueA}{RGB}{0,34,102}
\definecolor{blueB}{RGB}{0,43,128}
\definecolor{blueC}{RGB}{0,51,153}
\definecolor{blueD}{RGB}{0,60,179}
\definecolor{blueE}{RGB}{0,68,204}
\definecolor{blueF}{RGB}{0,77,230}




\definecolor{greyDDD}{rgb}{0.2,0.2,0.2}
\definecolor{greyDD}{rgb}{0.3,0.3,0.3}
\definecolor{greyD}{rgb}{0.4,0.4,0.4}
\definecolor{greyL}{rgb}{0.6,0.6,0.6}
\definecolor{greyLL}{rgb}{0.7,0.7,0.7}
\definecolor{greyLLL}{rgb}{0.8,0.8,0.8}

\definecolor{codegreen}{rgb}{0.0,0.258,0.0}
\definecolor{codeblue}{rgb}{0.0,0.199,0.6}
\definecolor{codepurple}{rgb}{0.398,0.0,1.0}




%%------------------------------------------------------------------------------
%%------------------------------------------------------------------------------
