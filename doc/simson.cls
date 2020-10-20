%$HeadURL$
%$LastChangedDate$
%$LastChangedBy$
%$LastChangedRevision$
\NeedsTeXFormat{LaTeX2e}[2001/01/01]
\ProvidesClass{simson}
              [2007/05/31 v1.0 $Revision$ Document class for Simson user guides]
\RequirePackage{graphicx}
\RequirePackage{ifthen}
\RequirePackage{color}
\RequirePackage[latin1]{inputenc}
%\RequirePackage{ifpdf}
\RequirePackage[swedish,english]{babel}
%%
\DeclareOption{a4paper}{\PassOptionsToClass{a4paper}{memoir}\afour}
\DeclareOption{10pt}{\PassOptionsToClass{10pt}{memoir}}
\DeclareOption{11pt}{\PassOptionsToClass{11pt}{memoir}}
\DeclareOption{12pt}{\PassOptionsToClass{12pt}{memoir}}
%%
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{memoir}}
%%
\newcommand{\afour}{
\renewcommand{\@kthunitlength}{1.0cm}
}
\newcommand{\@kthfrontfont}{\sffamily}
\newcommand{\@kthfrontfontbold}{\sffamily\bfseries}
\newcommand{\@kthchapfont}{\sffamily}
\newcommand{\@kthchapnumfont}{\sffamily}
\newcommand{\@kthunitlength}{1cm}
\newcommand{\@kthcoverpicturewidth}{21.0cm}
%%
\ExecuteOptions{a4paper}
%%
\ProcessOptions\relax
%%
\LoadClass[twoside,openright,onecolumn,10pt]{memoir}
%%
%% Give additional pdf information
%%
\AtBeginDocument{
\ifpdf
\pdfinfo
{
 /Producer (Simson LaTeX class $Revision$)
 /CreationDate (\today)
}
\fi
%%
%% Redefine lengths
%%
\setstocksize{297mm}{210mm}
\settrimmedsize{297mm}{210mm}{*}
\setlength{\trimtop}{0pt}
\setlength{\trimedge}{\stockwidth}
\addtolength{\trimedge}{-\paperwidth}
\addtolength{\uppermargin}{-10mm}
\addtolength{\lowermargin}{-10mm}
\addtolength{\textheight}{20mm}
\setlrmargins{*}{*}{1.618}
\setulmargins{*}{*}{*}
\setlength{\marginparwidth}{1.3cm}
\setlength{\marginparsep}{0.5cm}
\setlength{\parindent}{0pt}
\setlength{\parskip}{1ex plus 0.5ex minus 0.2ex}
\sideparswitchfalse
\checkandfixthelayout
\pagenumbering{roman}
\@pageheadings
\pagestyle{kthstyle}
%\chapterstyle{demo}
\chapterstyle{kthchapstyle}
} % End of AtBeginDocument
%%
\AtEndDocument{%
\relax
} % End of AtEndDocument
%%
\newcommand{\@pageheadings}{%
\makepagestyle{kthstyle}
  \makeevenhead{kthstyle}
    {\sffamily\@kthno}
    {\sffamily\@kthheadertitle}
    {\sffamily\@kthissue}
  \makeoddhead{kthstyle}
    {\sffamily\@kthissue}
    {\sffamily\@kthheadertitle}
    {\sffamily\@kthno}
  \makeevenfoot{kthstyle}{}{\sffamily\thepage}{}
  \makeoddfoot{kthstyle}{}{\sffamily\thepage}{}
%%
\makepagestyle{kthchapstyle}
  \makeevenhead{kthchapstyle}
    {}
    {}
    {}
  \makeoddhead{kthchapstyle}
    {}
    {}
    {}
  \makeevenfoot{kthchapstyle}{}{\sffamily\thepage}{}
  \makeoddfoot{kthchapstyle}{}{\sffamily\thepage}{}
%%
\aliaspagestyle{chapter}{kthchapstyle}
%\aliaspagestyle{cleared}{kthstyle}
\aliaspagestyle{cleared}{empty}
\makeheadrule{kthstyle}{\textwidth}{1pt}
\setsubsecheadstyle{\large\sffamily\raggedright}
\setsecheadstyle{\Large\sffamily\raggedright}
\setsubsubsecheadstyle{\normalsize\sffamily\raggedright}
\setparaheadstyle{\normalsize\sffamily}
\setparaheadstyle{\normalsize\sffamily}
\renewcommand{\cftchapterfont}{\sffamily}
\renewcommand{\cftsectionfont}{\sffamily}
\renewcommand{\cftsubsectionfont}{\sffamily}
\renewcommand{\cftsubsubsectionfont}{\sffamily}
\renewcommand{\cftparagraphfont}{\sffamily}
\renewcommand{\cftchapterpagefont}{\bfseries\sffamily}
\renewcommand{\cftsectionpagefont}{\sffamily}
\renewcommand{\cftsubsectionpagefont}{\sffamily}
\renewcommand{\cftsubsubsectionpagefont}{\sffamily}
\renewcommand{\cftparagraphpagefont}{\sffamily}
%\chapterstyle{section}
\captionnamefont{\small\sffamily}
\captiontitlefont{\small\sffamily}
%%
\newlength{\numberheight}
\newlength{\barlength}
\makechapterstyle{kthchapstyle}{%
\setlength{\beforechapskip}{40pt}
\setlength{\midchapskip}{0pt}
\setlength{\afterchapskip}{40pt}
\abstractnamefont{\flushleft\normalfont\HUGE\@kthchapfont}
\abstracttextfont{\normalfont\normalsize}
\renewcommand{\chapnamefont}{\normalfont\LARGE\flushright\@kthchapfont}
\renewcommand{\chapnumfont}{\normalfont\HUGE\@kthchapnumfont}
\renewcommand{\chaptitlefont}{\vspace{-16mm}\flushleft\normalfont\HUGE\@kthchapfont}
\renewcommand{\chaptername}{\relax}
\renewcommand{\appendixname}{\relax}
\renewcommand{\printchaptername}{%
\chapnamefont\MakeUppercase{\@chapapp}}
\renewcommand{\chapternamenum}{}
\setlength{\numberheight}{18mm}
\setlength{\barlength}{\paperwidth}
\addtolength{\barlength}{-\textwidth}
\addtolength{\barlength}{-\spinemargin}
\addtolength{\barlength}{20mm}
\renewcommand{\printchapternum}{%
\makebox[0pt][l]{%
\hspace{.8em}%
\resizebox{!}{\numberheight}{\chapnumfont \thechapter}%
\hspace{.8em}%
\color{\@kthchapterbarcolor}\rule{\barlength}{\numberheight}\color{black}
}
}
}}% End of pageheadings definition
%%
%% Macros
%%
\setverbatimfont{\scriptsize\ttfamily}
%\setlength{\verbatimindent}{3cm}
%%
\newcommand{\@kthcoverpicture}{}
\newcommand{\kthcoverpicture}[1]{\renewcommand{\@kthcoverpicture}{\includegraphics[width=\@kthcoverpicturewidth]{#1}}}
%%
\newcommand{\@kthno}{TRITA-MEK 2007:07}
\newcommand{\kthno}[1]{\renewcommand{\@kthno}{#1}}
%%
\newcommand{\@kthdate}{\@date}
\newcommand{\kthdate}[1]{\renewcommand{\@kthdate}{#1}}
%%
\newcommand{\@kthdivision}{KTH Mechanics}
\newcommand{\kthdivision}[1]{\renewcommand{\@kthdivision}{#1}}
\newcommand{\@kthdivisionswe}{Institutionen för Mekanik, KTH}
\newcommand{\kthdivisionswe}[1]{\renewcommand{\@kthdivisionswe}{#1}}
%%
\newcommand{\@kthzipcode}{SE-100 44 STOCKHOLM}
\newcommand{\kthzipcode}[1]{\renewcommand{\@kthzipcode}{#1}}
\newcommand{\@kthzipcodeswe}{SE-100 44 STOCKHOLM}
\newcommand{\kthzipcodeswe}[1]{\renewcommand{\@kthzipcodeswe}{#1}}
%%
\newcommand{\@kthphone}{Tel: + 46 8 790 6000}
\newcommand{\kthphone}[1]{\renewcommand{\@kthphone}{Tel: #1}}
%%
\newcommand{\@kthissue}{\mbox{}}
\newcommand{\kthissue}[1]{\renewcommand{\@kthissue}{#1}}
%%
\newcommand{\@kthheadertitle}{\mbox{}}
\newcommand{\kthheadertitle}[1]{\renewcommand{\@kthheadertitle}{#1}}
%%
\newcommand{\@kthinnerbackone}{\null}
\newcommand{\kthinnerbackone}[1]{\renewcommand{\@kthinnerbackone}{#1}}
%%
\newcommand{\@kthinnerbacktwo}{\null}
\newcommand{\kthinnerbacktwo}[1]{\renewcommand{\@kthinnerbacktwo}{#1}}
%%
\newcommand{\@kthisbn}{\null}
\newcommand{\kthisbn}[1]{\renewcommand{\@kthisbn}{#1}}
%%
\newcommand{\@kthchapterbarcolor}{red}
\newcommand{\kthchapterbarcolor}[1]{\renewcommand{\@kthchapterbarcolor}{#1}}
%%
\newcommand{\marginlabel}[1]%
{\sidepar{\raggedright{\sffamily\scriptsize #1}}}
%%
%% Frontpage definition
%%
\renewcommand\maketitle{%
  \pagestyle{empty}\relax%
  \frontpage
  \cleartooddpage%
  \setcounter{page}{1}
  \innerpage
  \clearpage
  \innerbackpage
  \clearpage
  \pagestyle{kthstyle}\relax%
}
%%
%% Front page definition
%%
\newcommand{\frontpage}{%
\thispagestyle{empty}%
\newsavebox{\descriptionfrontbox}%
\setlength{\unitlength}{\@kthunitlength}
\sbox{\descriptionfrontbox}{%
%%
%% Front page layout
%%
\setlength{\unitlength}{1.0cm}
\begin{picture}(19.,26.)
%%
%% KTH logo
%%
\put(-0.9,23.5){\includegraphics[width=30 mm]{figures/KTH_eng_CMYK}}
%%
%% Title
%%
\put(4.6,27.3){\parbox[t]{13.8cm}{\begin{flushright}\HUGE\sffamily\@title\end{flushright}}}
%%
%% Author
%%
\put(-0.85,21.3){\parbox[t]{10.2cm}{\sffamily\LARGE\@author}}
%%
%% Coverpicture
%%
\put(-1.40,3.7){\@kthcoverpicture}
%%
%% KTH Mechanics text
%%
\put(-0.85,0.7){\makebox[8.0cm][l]{\sffamily\bfseries{\@kthdivision}}}
\put(-0.85,0.2){\makebox[8.0cm][l]{\sffamily\@kthno}}
%%
\end{picture}}%
\null\hspace{-1cm}
\begin{picture}(0,0)
\put(-0.9,-25.0){\usebox{\descriptionfrontbox}}
\end{picture}}
%%
%% Inner front page layout
%%
\newcommand{\innerpage}{%
\thispagestyle{empty}%
\newsavebox{\descriptioninnerbox}%
\setlength{\unitlength}{\@kthunitlength}
\sbox{\descriptioninnerbox}{%
%%
%% Inner page layout
%%
\setlength{\unitlength}{1.0cm}
\begin{picture}(19.,26.)
%%
%% Title
%%
\put(1.6,24.0){\parbox[t]{14.0cm}{\begin{center}\HUGE\sffamily\@title\end{center}}}
%%
%% Author
%%
\put(1.6,19.0){\parbox[t]{14.0cm}{\begin{center}\sffamily\Huge\@author\end{center}}}
%%
%% KTH Mechanics text
%%
\put(-0.85,0.7){\makebox[8.0cm][l]{\sffamily\bfseries{\@kthdivision}}}
\put(-0.85,0.2){\makebox[8.0cm][l]{\sffamily\@kthno}}
%%
\end{picture}}%
\null\hspace{-1cm}
\begin{picture}(0,0)
\put(-0.9,-25.0){\usebox{\descriptioninnerbox}}
\end{picture}}
%%
%% Back side of inner front page
%%
\newcommand{\innerbackpage}{%
\thispagestyle{empty}%
\newsavebox{\descriptioninnerbackbox}%
\setlength{\unitlength}{\@kthunitlength}
\sbox{\descriptioninnerbackbox}{%
%%
%% Inner back page layout
%%
\setlength{\unitlength}{1.0cm}
\begin{picture}(19.,26.)
%%
%% Copyright statement
%%
\put(0.0,23.5){\parbox[t]{13.0cm}{This is software which is
distributed freely on a limited basis; it comes with no guarantees
whatsoever. Problems can be reported to henning@mech.kth.se, but no
action is promised.  If results obtained by using the programs
included in the \esoft{Simson} distribution are published the authors
would like an acknowledgment, \emph{e.g.}\ in the form of referencing 
this code manual.}}
%%
%% Front image description
%%
\put(0.0,21.5){\parbox[t]{13.0cm}{\@kthinnerbackone}}
%%
%% Revision information
%%
\put(0.0,18.0){\parbox[t]{13.0cm}{\@kthinnerbacktwo}}
%%
%% ISBN
%%
\put(0.0,17.2){\parbox[t]{13.0cm}{\@kthisbn}}
%%
\end{picture}}%
\null\hspace{-1cm}
\begin{picture}(0,0)
\put(-0.9,-25.0){\usebox{\descriptioninnerbackbox}}
\end{picture}}
%%
%% Set toc display depth
%%
\maxsecnumdepth{subsubsection}
\setsecnumdepth{subsubsection}
\maxtocdepth{subsection}
\settocdepth{subsection}
%%
\definecolor{green}{cmyk}{0.18,0.0,1.0,0.31}
\definecolor{red}{cmyk}{0.0,1.0,0.64,0.33}
%%\definecolor{blue}{cmyk}{1.0,0.83,0.46,0.13}
\definecolor{lilac}{cmyk}{0.65,0.86,0.49,0.0}
\definecolor{lightblue}{cmyk}{0.85,0.24,0.0,0.0}
\definecolor{gray}{cmyk}{0.0,1.0,0.0,0.51}
\definecolor{black}{cmyk}{0.0,0.0,0.0,1.0}
\definecolor{white}{cmyk}{0.0,0.0,0.0,0.0}
%%
%% To make a correct alphabetic order of authors in the bibliography also
%% when having names including the Swedish letters åäö
%%
\newcommand{\noopsort}[1]{}
%%
%\newcommand{\td}{{\text d\rm}}
%\newcommand\etal{\mbox{\emph{et al.}}}
%\newcommand\etc{etc.\ }
%\newcommand\eg{e.g.\ }
%%
%% Keyword tags
%%
\newcommand{\epar}[1]{\textbf{#1}}                                       % Parameter text
\newcommand{\esoft}[1]{\textbf{{#1}}}                                    % Software font
\newcommand{\efile}[1]{\textsf{\textbf{{\small #1\normalsize}}}}         % File font
%\newcommand{\efile}[1]{\textsf{\textbf{{\normalsize #1\normalsize}}}}    % File font
\newcommand{\efilehead}[1]{\textsf{\textbf{#1}}}                         % File section font
\newcommand{\epath}[1]{\textsf{\textbf{{\small #1\normalsize}}}}         % Path
\newcommand{\esub}[1]{\textbf{#1}}                                       % Subroutine
\newcommand{\envvar}[1]{\textsf{{\small #1\normalsize}}}                 % Environment variable
\newcommand{\eopt}[1]{\texttt{{\small #1\normalsize}}}                   % Program options
\newcommand{\ecommand}[1]{\texttt{{#1}}}                                 % Command
%%
\setlength{\parskip}{0pt}
\endinput
