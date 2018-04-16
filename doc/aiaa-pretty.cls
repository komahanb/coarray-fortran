% I dedicate this template unto thee, Sara C. Spangelous.
%
%
%%
%% This is file `aiaa-pretty.cls'.  It is based off of aiaa-tc.cls, which was
%% written by Bil Kleb, Bill Wood, and Erich Knauseberger.
%%
%% Editing this file is perfectly permitted by the distributing author, Derek
%% Dalle, but changing the name is recommended before making any major changes.
%%
%% -----------------------------------------------------------------------------
%% IMPORTANT USE GUIDELINES:
%%
%% This document class works like the standard 'article' class except for the
%% way it handles abstracts.  Rather than putting the text of the abstract in an
%% environment, the user now specifies the text of the abstract as the argument
%% of a command.
%%
%% Correct:
%%
%%   \abstract{This is a paper about ...}
%%
%% Incorrect:
%%
%%   \begin{abstract}
%%      This is a paper about ...
%%   \end{abstract}
%%
%% -----------------------------------------------------------------------------
%% This provides a template for AIAA papers.
%% 
%% This software is released under the terms of the LaTeX Project Public
%% License.  Copyright (C) 2004 by Bil Kleb, Bill Wood, and Erich Knausenberger.
%%
%% Description: LaTeX package for AIAA papers.
%% Keywords: LaTeX, class, AIAA, BibTeX, bibliographic-style
%% Author (aiaa-pretty.cls): Derek Dalle
%% Author (aiaa-tc.cls): Bil Kleb, Bill Wood, Erich Knausenberger
%% Maintainer: Derek Dalle
%% Version: 2.0 <25 May 2011>
%%
%% There are two installation options.
%%
%%   1. Put this file in the same directory as your .tex file, or near it, and 
%%      specify the relative path to this file from the .tex folder in the 
%%      \documentclass line.  For example, if this file is in a directory one
%%      level higher than your .tex file, you should use
%%      \documentclass{../aiaa-pretty}.
%%
%%   2. Alternately, you could put this file somewhere in your TeX distribution
%%      directory.  If 'texmf' is the path to your main TeX distribution folder,
%%      then a good recommendation is putting this file in 
%%      'texmf/tex/latex/aiaa/'.  This is the same directory that would contain
%%      the aiaa-tc.cls file if it is installed.
%%
%%  'aiaa-pretty.cls' - A LaTeX class to create AIAA technical papers.
%%
%%  The basic approach is to load the LaTeX article class and various
%%  standard LaTeX packages, and then make modifications and extensions to
%%  conform to AIAA's guidelines.
%%
%%  To use, simply select this class file in your LaTeX document.  For example,
%%
%%    \documentclass{aiaa-pretty}
%%
%%  Class options:
%%
%%   backrefs   - using this option puts PDF links in the bibliography that
%%                take the reader to the page(s) on which the reference is
%%                used.
%%
%%   conference - used for hand carrying a paper to a 'loose papers'
%%                conference.  It adds the conference name to the
%%                top of the first page, paper number to the bottom
%%                of each page, and a copyright notice on the bottom
%%                of the first.
%%
%%   journal    - produces a style that attempts to replicate a typical
%%                style found in an AIAA journal
%%
%%   plain      - this option removes any words or labels referring to to the
%%                AIAA.  It is most convenient when you are preparing either
%%                a draft or a paper for another agency.  You can also change
%%                the label to something else using the command
%%                     \footerlabel{NOT the AIAA}
%%
%%   submit     - produces output typically required for journal submission,
%%                i.e., 12pt font and 'double spacing'.
%%
%%  Additional options specified in the \documentclass line will be passed 
%%  to the article class.
%%
%%  Class-specific commands:
%%
%%   \thanksibid{} - Reuses \thanks footnote symbol from previous
%%                   author specified by position in \author command.
%%                   For example, `someone\thanksibid{2}' will reuse the
%%                   footnote symbol from the second author for 'someone'
%%                   without producing a repeated footnote at the bottom.
%%
%%   \AIAAauthor[short info]{full names, etc.}
%%                      - New method for specifying author information.
%%                        The old \author{names, etc.} is still functional,
%%                        but the headers for the 'journal' option will not
%%                        contain the authors' last names.
%%
%%   \AIAAabstract{text of abstract}
%%                      - New method of entering abstract.  The old
%%                        \begin{abstract}... text ... \end{abstract} method
%%                        is still functional, but the abstract will be in
%%                        one column when using the 'journal' option, which
%%                        match the AIAA journal formatting.  Make sure to
%%                        use this command BEFORE the \maketitle line.
%%
%%   \AIAApapernumber{} - Inserts number of AIAA paper.  For use with the
%%                        'conference' option.
%%
%%    \AIAAconference{} - Inserts the name of the conference.  For use with
%%                        the 'conference' option.
%%
%%       \AIAAjournal{} - Inserts the name of the journal.  For use with the
%%                        'journal' option.
%%
%%     \AIAAcopyright{} - Inserts copyright info on bottom of first page
%%                        when used with 'conference' option.  Must be filled
%%                        with one of the various copyright types:
%%                          \AIAAcopyrightA{year},
%%                          \AIAAcopyrightB{year}{holder},
%%                          \AIAAcopyrightC,
%%                          \AIAAcopyrightD{year}, or
%%                          \AIAAcopyrightDAPL{year}.
%%
%%   \AIAAsubmitinfo{} - Footer information used by 'submit' option.
%%
%%

%% ----- INITIAL DEFINITIONS ---------------------------------------------------
% This prevents the compiler from running on old versions of LaTeX.
\NeedsTeXFormat{LaTeX2e}

% This command gives the name of the class.
\ProvidesClass{aiaa-pretty}[2010/01/10 v2.0 AIAA paper]


%% ----- OPTIONS ---------------------------------------------------------------
% LaTeX help: logical variables
% Create a new logical variable using
%   \newif\if@something
% This creates a variable \@something with a default of false.

% These variables control the document type.
\newif\if@aiaa@conference
\newif\if@aiaa@submit
\newif\if@aiaa@journal
\newif\if@aiaa@plain

% These variables control which items to show.
\newif\if@q@abstract
\newif\if@q@author
\newif\if@q@shortauthor
\newif\if@q@title
\newif\if@q@footer
\newif\if@q@header
\newif\if@q@backref
\newif\if@q@copyright
\newif\if@q@papernumber
\newif\if@q@conference
\newif\if@q@lastpage
\newif\if@q@lastpage@link

% This value can be used to hide any AIAA labels.
\newif\if@q@plain

% These variables relate to the hyperref links.
\newif\if@q@links
\newif\if@links@plain
\newif\if@links@bright

% This option allows the user to compress the entire document.
\newif\if@q@compress

% This sets up any of the above variables that should default to true.
\@q@linkstrue
\@q@lastpagetrue

% Document format options
\DeclareOption{conference}{\@aiaa@conferencetrue}
\DeclareOption{submit}{\@aiaa@submittrue}
\DeclareOption{journal}{\@aiaa@journaltrue\@q@lastpagefalse}
% Option to hide footers
\DeclareOption{plain}{\@q@plaintrue}
% Option to compress document
\DeclareOption{compress}{\@q@compresstrue}
% Option to use backward references for references
\DeclareOption{backrefs}{\@q@backreftrue}
% Option to put the number of pages on each page.
\DeclareOption{lastpage}{\@q@lastpagetrue}
% Option to not put the number of pages on each page.
\DeclareOption{lastpage=false}{\@q@lastpagefalse}
% Option to make the number of pages a link to the last page.
\DeclareOption{lastpage=link}{\@q@lastpage@linktrue}
% Options to change the colors of the links
\DeclareOption{links=false}{\@q@linksfalse}
\DeclareOption{links=plain}{\@links@plaintrue}
\DeclareOption{links=bright}{\@links@brighttrue}

% This passes any other options on to the article class.
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}%
 \typeout{Unknown option: passing '\CurrentOption'
  option on to the standard LaTeX article class.}}

% This command appears to tell LaTeX to simply process the class options.
\ProcessOptions


%% ----- DOCUMENT FORMAT -------------------------------------------------------
% This control structure calls the basic LaTeX article.cls class.
\if@aiaa@submit
 % Submit class: large font and spacing in single column
 \LoadClass[letterpaper,10pt]{article}
 % Make sure that other formats are suppressed.
 \@aiaa@conferencefalse \@aiaa@journalfalse \@aiaa@plainfalse
\else\if@aiaa@journal
 % Journal class: mirror the appearance of an AIAA journal article
 \LoadClass[letterpaper,10pt,twoside]{article}
\else
 % Conference class: format appropriate for AIAA conference paper
 \LoadClass[letterpaper,10pt]{article}
\fi\fi

% This class uses the geometry.sty package to set the margins.
\if@aiaa@journal
 % The margins are very small for the journal version.
 \RequirePackage[margin=1.5cm,bottom=1.1cm,footskip=8pt,headsep=0.4cm,%
  headheight=\baselineskip]{geometry}
 % This sets the distance between columns.
 \setlength{\columnsep}{0.75cm}
\else
 % Larger margins for conference class
 \RequirePackage[margin=1in,top=0.75in,bottom=0.9in, %
  headheight=0.5\baselineskip]{geometry}
\fi


%% ----- PACKGAGES -------------------------------------------------------------
% This section tells LaTeX to call a bunch of other packages.
\RequirePackage{graphicx}   % images
\RequirePackage{array}      % extra row height
\RequirePackage{lastpage}   % to get reference for last page number
\RequirePackage{fancyhdr}   % simple control of headers and footers
\RequirePackage{multicol}   % package for multi-column nomenclature
\RequirePackage{mathptmx}   % changes font to Times New Roman
\RequirePackage{verbatim}   % New verbatim format
\RequirePackage{textcomp}   % proper single quotes in verbatim
\RequirePackage{ifthen}     % logical structure
% \RequirePackage[BF,bf]{subfigure}     % subfigures
\RequirePackage{amsmath}    % AMS packages
\RequirePackage{amsfonts}
\RequirePackage{amssymb}

% Need a package to set 1.5-spaced lines for 'submit'
\if@aiaa@submit\RequirePackage{setspace}\fi
% Need a package to fix the author footnotes in two-column format
\if@aiaa@journal\RequirePackage{abstract}\fi


%% ----- IDENTIFIERS AND ABSTRACT ----------------------------------------------
% Variables that hold the value of various identifiers.
\renewcommand{\@author}{}
\newcommand{\@shortauthor}{}
\newcommand{\@abstract}{}
\renewcommand{\@title}{}
\newcommand{\@copyright}{}
\newcommand{\@submitinfo}{}
\newcommand{\@conference}{}
\newcommand{\@journal}{}

% New author command
\renewcommand*{\author}[2][]{ %
 % Set the main author field.
 \renewcommand{\@author}{#2} %
 % Set a flag to use it.
 \@q@authortrue %
 % Check for an optional argument.
 \ifthenelse{\equal{#1}{}}{}{ %
  % Set the short author field.
  \renewcommand{\@shortauthor}{#1} %
  % Set a flag to use it.
  \@q@shortauthortrue %
 } %
}
% Command to set just the short author.
\newcommand{\shortauthor}[1]{ %
 \renewcommand{\@shortauthor}{#1} %
 \@q@shortauthortrue}
% Legacy command
\let\AIAAauthor\author

% New abstract command
\renewcommand{\abstract}[1]{ %
 \renewcommand{\@abstract}{} %
 \@q@abstracttrue}
% Legacy command
\let\AIAAabstract\abstract

% New title command
\renewcommand{\title}[1]{ %
 \renewcommand{\@title}{#1} %
 \@q@titletrue}
% Legacy command
\let\AIAAtitle\title

% New abstract command
\renewcommand{\abstract}[1]{
 \renewcommand{\@abstract}{#1}\@q@abstracttrue}
% Legacy command
\let\AIAAabstract\abstract

% These are counters for the author affiliations.
\newcounter{aiaa@cthanks}
\newcounter{aiaa@cfootnote}
% This command enables multiple authors to have the same footnote hook.
\newcommand{\thanksibid}[1]{%
\begingroup
% Store the current value of the footnote counter.
\setcounter{aiaa@cfootnote}{\value{footnote}}%
% Set the thanks counter to the argument
\setcounter{aiaa@cthanks}{#1}%
% Get the footnote symbol corresponding to above
\renewcommand\thefootnote{\fnsymbol{aiaa@cthanks}}%
% Make a footnote.
\footnotemark%
% Reset the footnote counter to the original value.
\setcounter{footnote}{\value{aiaa@cfootnote}}%
\endgroup %
}

% Commands to insert the identifiers
\newcommand{\insertauthor}{\@author}
\newcommand{\insertshortauthor}{\@shortauthor}
\newcommand{\inserttitle}{\@title}
\newcommand{\insertabstract}{\@abstract}


%% ----- FLOATS ----------------------------------------------------------------
% This eliminates the parentheses around the sublabel in references.
% \renewcommand{\thesubfigure}{\alph{subfigure}}
% This formats the actual label area for each subfigure.
% \renewcommand{\@thesubfigure}{\alph{subfigure})\hskip\subfiglabelskip}
% This eliminates the parentheses around the sublabel in references.
% \renewcommand{\thesubtable}{\alph{subtable}}
% This formats the actual label area for each subtable.
% \renewcommand{\@thesubtable}{\alph{subtable})\hskip\subfiglabelskip}

% This sets the distance between a float and the paragraph.
\setlength{\intextsep}{0.7\baselineskip \@plus 2\p@ \@minus 1\p@}

% This sets the priorization for where floats should go.
\renewcommand{\fps@figure}{htbp}
\renewcommand{\fps@table}{htbp}

% Shorten the length below the caption
\setlength\belowcaptionskip{5pt plus 1pt minus 1pt}
% Change the figure caption styles.
\renewcommand{\@makecaption}[2]{%
 % Space between caption and figure
 \vskip\abovecaptionskip %
 % Actual caption insertion
 \sbox\@tempboxa{\footnotesize\bfseries #1.~~#2} %
 \ifdim \wd\@tempboxa >\hsize
  \footnotesize\bfseries #1.~~#2\par
 \else
  \global \@minipagefalse
  \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
 \fi
 % Space between caption and next object
 \vskip\belowcaptionskip}
 

%% ----- LENGTHS AND SIZES -----------------------------------------------------
% Vertical spaces for the title area
\newlength{\toptitleskip}
\newlength{\titleauthorskip}
\newlength{\titleabstractskip}
\newlength{\authorabstractskip}
\newlength{\titletextskip}

% Set the lengths of the title area skips
\setlength{\toptitleskip}{1ex}
\setlength{\titleauthorskip}{5ex}
\setlength{\titleabstractskip}{5ex}
\setlength{\authorabstractskip}{5ex}
\setlength{\titletextskip}{5ex plus 2pt minus 3pt}

% Command to make the skips in the title area smaller.
\newcommand{\compresstitle}{ %
 \setlength{\toptitleskip}{-1ex plus 2ex}
 \setlength{\titleauthorskip}{2ex}
 \setlength{\titleabstractskip}{2ex}
 \setlength{\authorabstractskip}{2ex}
 \setlength{\titletextskip}{3ex plus 2pt minus 3pt}
}

% Command to make the section skips smaller.
\newcommand{\compresssections}{ %
 % This gives formatting rules for section headers.
 \renewcommand\section{%
  \@startsection {section}{1}{\z@}%
   {0.1ex \@plus -1ex \@minus -.2ex}%
   {0.7ex \@plus.2ex}%
   {\center\normalfont\large\bfseries}}
 % This gives formatting rules for subsection headers.
 \renewcommand\subsection{%
  \@startsection{subsection}{2}{\z@}%
   {-5.25ex\@plus -1ex \@minus -.2ex}%
   {0.5ex \@plus .2ex}%
   {\normalfont\normalsize\bfseries}}
}

% Command to set the font size to 9 with normal spacing
% It is necessary because there is no 9pt option for article.
\newcommand{\journalfont}{ %
 % This does the main font resizing work.
 \renewcommand\normalsize{ %
  % The first number is the size of the text.
  % The second number is the size of the line.
  \@setfontsize\normalsize\@ixpt\@xpt%
  % Set the amount of vertical white space before and after equations.
  \abovedisplayskip 8\p@ \@plus4\p@ \@minus4\p@
  \belowdisplayskip 8\p@ \@plus4\p@ \@minus4\p@
  \abovedisplayshortskip \z@ \@plus3\p@
  \belowdisplayshortskip 5\p@ \@plus3\p@ \@minus3\p@
  % I do not know what this means.
  \let\@listi\@listI
 }
 % The small size needs tweaking, too
 \renewcommand\small{%
  % Decrease the font size to 8pt.
  \@setfontsize\small\@viiipt{9.5}%
  % Decrease the margins above and below equations.
  \abovedisplayskip 6\p@ \@plus2\p@ \@minus4\p@
  \abovedisplayshortskip \z@ \@plus\p@
  \belowdisplayshortskip 3\p@ \@plus\p@ \@minus2\p@
  % I still don't know what this means.
  \def\@listi{\leftmargin\leftmargini
   \topsep 3\p@ \@plus\p@ \@minus\p@
   \parsep 2\p@ \@plus\p@ \@minus\p@
   \itemsep \parsep} %
 }
}

% Command to set the font size to 9 with tight spacing
% It is necessary because there is no 9pt option for article.
\newcommand{\smallfonts}{ %
 % This does the main font resizing work.
 \renewcommand\normalsize{ %
  % The first number is the size of the text.
  % The second number is the size of the line.
  \@setfontsize\normalsize\@ixpt\@ixpt%
  % Set the amount of vertical white space before and after equations.
  \abovedisplayskip 4\p@ \@plus4\p@ \@minus2\p@
  \belowdisplayskip 4\p@ \@plus4\p@ \@minus2\p@
  \abovedisplayshortskip \z@ \@plus3\p@
  \belowdisplayshortskip 3\p@ \@plus3\p@ \@minus2\p@
  % I do not know what this means.
  \let\@listi\@listI
 }
 % The small size needs tweaking, too
 \renewcommand\small{%
  % Decrease the font size to 8pt.
  \@setfontsize\small\@viiipt{9}%
  % Decrease the margins above and below equations.
  \abovedisplayskip 3\p@ \@plus2\p@ \@minus1\p@
  \abovedisplayshortskip \z@ \@plus\p@
  \belowdisplayshortskip 2\p@ \@plus\p@ \@minus\p@
  % I still don't know what this means.
  \def\@listi{\leftmargin\leftmargini
   \topsep 3\p@ \@plus\p@ \@minus\p@
   \parsep 2\p@ \@plus\p@ \@minus\p@
   \itemsep \parsep} %
 }
}


%% ----- TITLE AREA ------------------------------------------------------------
\newcommand{\tskrp}{5ex}
%\setlength\titleauthorskip{12ex plus 2pt minus 3pt}

% Define a new command that does most of the work of making the title.
\renewcommand{\@maketitle}{ %
 % Housekeeping
 \newpage %
 % Change the footnote style to silly marks.
 \renewcommand\thefootnote{\@fnsymbol\c@footnote} %
 % Insert the information in a centered environment
 \begin{center} %
  % Initial vertical space
  \vskip\toptitleskip %
  % Set the footnote command equal to the thanks command.
  \let\footnote\thanks %
  % Check for a title
  \if@q@title %
   {\LARGE \bfseries \inserttitle} %
  \fi %
  % Check for an author field.
  \if@q@author %
   % Title-to-author gap
   \if@q@title \vskip\titleauthorskip \fi %
   % Insert the author.
   { %
    % Increase font size for all fields of author.
    \large\sloppy %
    % Insert author in a tabular environment.
    \hfill\begin{tabular}[t]{c}\insertauthor\end{tabular}\hfill
   } %
  \fi %
  % Check for an abstract.
  \if@q@abstract %
   % Check which gap to insert above.
   \if@q@author %
    % Author-to-abstract skip
    \vskip\authorabstractskip %
   \else\if@q@title %
    % Title-to-abstract skip
    \vskip\titleabstractkip %
   \fi\fi %
   % Insert the abstract
   \begin{minipage}{13.5cm} %
    \bfseries \small \insertabstract %
   \end{minipage} %
   % Abstract-to-text gap
   \vskip\titletextskip %
  \else\if@q@author %
   % Author-to-text gap
   \vskip\titletextskip %
  \else\if@q@title %
   %Title-to-text gap
   \vskip\titletextskip %
  \fi\fi\fi %
 \end{center} %
 % Produce the author footnotes.
 \@thanks %
 % Restart the footnotes (Don't count the 'thanks' information.).
 \setcounter{footnote}{0} %
 % Set the footnote symbol to letters.
 \renewcommand\thefootnote{\@alph\c@footnote} %
}

% This section distinguishes between the two-column and one-column options.
\if@aiaa@journal %
 % For a journal, use two-column format except for the title and abstract.
 \renewcommand{\maketitle}{ %
  % Two-column title area
  \twocolumn[\@maketitle] %
  % Remove extra headers and footers from titlepage.
  \thispagestyle{plain} %
  % This is needed to print the footnotes.
  \saythanks %
  % Reset the footnote style.
  \renewcommand\thefootnote{\@alph\c@footnote}}
\else
 % Insert the title area as part of the whole one-column document.
 \renewcommand{\maketitle}{\@maketitle}
\fi


%% ----- HEADER AND FOOTER -----------------------------------------------------
% Whether or not there is a footer label to show
\newif\if@q@footer@label
% Commands to hide or show footer label
\newcommand{\showfooterlabel}{\@q@footer@labeltrue}
\newcommand{\hidefooterlabel}{\@q@footer@labelfalse}
% Command to hold the value of the footer
\newcommand{\@aiaa@footer@label}{ %
 American Institute of Aeronautics and Astronautics}
% Command to set it
\newcommand{\footerlabel}[1]{ %
 \renewcommand{\@aiaa@footer@label}{#1} \showfooterlabel}
% Command to insert it
\newcommand{\insertfooterlabel}{\@aiaa@footer@label}

% Set the default value of the footer label.
\if@q@plain
\else
 \footerlabel{American Institute of Aeronautics and Astronautics}
\fi

% Command to hold the value of the page reference
\newcommand{\@pagenum}{}
% Command to reset the page reference value
\newcommand{\@refresh@pagenum}{ %
 \if@q@lastpage
  % Check if links are turned on.
  \if@q@links
   % Check if the option for last-page links is on.
   \if@q@lastpage@link
    % Link to last page
    \renewcommand{\@pagenum}{\thepage~of \pageref{LastPage}}
   \else
    % Last page but no link
    \renewcommand{\@pagenum}{\thepage~of \pageref*{LastPage}}
   \fi
  \else
   % Last page but no links in document
   \renewcommand{\@pagenum}{\thepage~of \pageref{LastPage}}
  \fi
 \else
  % Nothing at all about last page
  \renewcommand{\@pagenum}{\thepage}
 \fi}
% Run this command
\@refresh@pagenum
% Command to get rid of last page.
\newcommand{\hidelastpage}{\@q@lastpagefalse\@refresh@pagenum}
% Command to show last page.
\newcommand{\showlastpage}{\@q@lastpagetrue\@refresh@pagenum}
% Command to show last page and link it.
\newcommand{\linklastpage}{ %
 \@q@lastpagetrue\@q@lastpage@linktrue %
 \@refresh@pagenum}


% This section sets up the footer.
\if@aiaa@journal
 % This eliminates the line at the top of pages.
 \renewcommand{\headrulewidth}{0pt}
 % This redefines the "plain" page style, which is used on the first page.
 \fancypagestyle{plain}{ %
  % First clear out the regular page style
  \fancyhf{} %
  % Put the title of the journal in the upper left.
  \fancyhead[L]{\scshape\@journal} %
  % Put the page number (1) in the bottom center
  \fancyfoot[C]{\thepage}} %
 % Clear the page style for all other pages.
 \fancyhf{} %
 % Page number on upper outer corner
 \fancyhead[RE,LO]{\footnotesize\@pagenum} %
 % Authors' last names on top center
 \fancyhead[CE,CO]{\footnotesize\MakeUppercase{\insertshortauthor}}
\else
 % This clears out the regular page style.
 \fancyhf{} %
 % This turns off the line at the top of pages.
 \renewcommand{\headrulewidth}{0pt} %
 % This sets the value of the footer for each page.
 \cfoot{
  % Set the font size and type.
  \footnotesize\normalfont
  % Check for a last page reference.
  \@pagenum \\
  % Check if there is a footer label.
  \if@q@footer@label %
   % Put a small vertical under the page number.
   \rule[.2\baselineskip]{0.5in}{0.2pt} \\ %
   % Insert the footer label.
   \insertfooterlabel %
  \fi } %
\fi
% This tells LaTeX to actually use the customized page style.
\pagestyle{fancy}


%% ----- SECTION HEADERS -------------------------------------------------------
% Update section numbering.
\renewcommand\thesection      {\@Roman\c@section}
\renewcommand\thesubsection   {\@Alph\c@subsection}
\renewcommand\thesubsubsection{\@arabic\c@subsubsection}
\renewcommand\theparagraph    {\@alph\c@paragraph}
\renewcommand\thesubparagraph {\@roman\c@subparagraph}
% No counters \paragraph and below
\setcounter{secnumdepth}{3}
% This confusing command puts a period after the section numbers.
\renewcommand{\@seccntformat}[1]{\csname the#1\endcsname.\quad}

% This gives formatting rules for section headers.
\renewcommand\section{%
 \@startsection {section}{1}{\z@}%
  {3.5ex \@plus 1ex \@minus -.2ex}%
  {1.0ex \@plus.2ex}%
  {\center\normalfont\large\bfseries}}

% This gives formatting rules for subsection headers.
\renewcommand\subsection{%
 \@startsection{subsection}{2}{\z@}%
  {-3.25ex\@plus -1ex \@minus -.2ex}%
  {0.5ex \@plus .2ex}%
  {\normalfont\normalsize\bfseries}}

% This gives formatting rules for subsubsection headers.
\renewcommand\subsubsection{%
 \@startsection{subsubsection}{3}{\z@}%
  {-3.25ex\@plus -1ex \@minus -.2ex}%
  {1.5ex \@plus .2ex}%
  {\normalfont\normalsize\itshape}}

% If you are silly enough to use paragraphs with headings, this controls them.
\renewcommand\paragraph{%
 \@startsection{paragraph}{4}{\z@}%
  {2.25ex \@plus1ex \@minus.2ex}%
  {-1em}%
  {\normalfont\normalsize\scshape}}

% Formatting for subparagraphs.
\renewcommand\subparagraph{%
 \@startsection{subparagraph}{5}{\parindent}%
  {2.25ex \@plus1ex \@minus .2ex}%
  {-1em}%
  {\normalfont\normalfont}}


%% ----- APPENDIX --------------------------------------------------------------
% Actual coding for the AIAA-style appendices
\newcommand{\aiaaappendix}[1]{\section*{Appendix: #1}}


%% ----- BIBILIOGRAPHY ---------------------------------------------------------
% Save the original bibliography command.
\let\@tex@bibliography\bibliography

% Create a new command for the bibliography.
\renewcommand*{\bibliography}[1]{ %
 % Update the page format.
 \if@aiaa@submit
  % Switch to single spacing
  \singlespacing %
 \fi %
 % Use a small font size.
 \small %
 % Create a pdf anchor.
 \phantomsection %
 % Add to table of contents.
 \addcontentsline{toc}{chapter}{References} %
 % Load the bibliography.
 \@tex@bibliography{#1} %
}


%% ----- CODE ENVIRONMENT ------------------------------------------------------
% This generates a 'code' environment that allows inline citations
% of small snippets of code.
\newenvironment{code} %
 {\quote\small\verbatim} %
 {\endverbatim\endquote}

% Formatting commands
% This command is used for names of files
\newcommand{\tfile}[1]{{\footnotesize\textsf{#1}}}
% This command puts a box around text that is meant for menu entries.
\newcommand{\tmenu}[1]{\framebox[\width]{\ {\tiny\textsf{#1}}\ }}
% This command is intended for names of variables within programs.
\newcommand{\tvar}[1]{{\footnotesize\texttt{#1}}}
% This command is used to display the contents of strings.
% It automatically adds the correct, upright quotes around the argument.
\newcommand{\tstring}[1]{ %
 {\footnotesize\texttt{\textquotesingle#1\textquotesingle}}}
% This command gives a bold version of the variable names.
\newcommand{\tfunction}[1]{{\footnotesize\texttt{\textbf{#1}}}}

%% ----- MATH COMMANDS ---------------------------------------------------------
% Command for Leibniz notation of derivatives
\newcommand{\dd}[2]{\frac{d#1}{d#2}}
% Command for Leibniz notation of partial derivatives
\newcommand{\pd}[2]{\frac{\partial#1}{\partial#2}}
% Command for norm
\newcommand{\norm}[1]{\lVert#1\rVert}
% Command for absolute value
\newcommand{\abs}[1]{\lvert#1\rvert}


%% ----- TEXT COMMANDS ---------------------------------------------------------
% Command for 'and others'
\newcommand{\etal}{\textit{et al.}}
% Command for 'for example'
\newcommand{\eg}{\textit{e.g.}}
% Command for 'that is'
\newcommand{\ie}{\textit{i.e.}}
% Command for 'and the rest'
\newcommand{\etc}{\textit{etc.}}

% Short space
\newcommand{\sspace}{\hspace{0.18em}}

%% ----- NOMENCLATURE ENVIRONMENT ----------------------------------------------
% New command for nomenclature environment
\newenvironment{nomenclature}[1][]
{ %
 % Insert the title if present.
 \ifthenelse{\equal{#1}{}}{ %
  % Start the tabular environment.
  \noindent\tabular{@{} p{3.5em} @{\hspace{2pt}=\hspace{1em}} l} %
 }{ %
 % Insert the title.
 \noindent\textit{#1}\\[3pt] %
 % Start the tabular environment.
 \noindent\tabular{@{} p{3.5em} @{\hspace{2pt}=\hspace{1em}} l} %
 }
}{ %
 % End the tabular environment
 \endtabular %
}

%% ----- LINKS -----------------------------------------------------------------
% Check if links should be loaded at all
\if@q@links
 % Load the link package.
 \if@q@backref
  % Load the package with backwards references for bibliography.
  \RequirePackage[pagebackref=true,colorlinks=true]{hyperref}
 \else
  % No back references
  \RequirePackage[colorlinks=true]{hyperref}
 \fi
 % Check the colors of the links that should be used.
 \if@links@plain
  % Use all black links.
  \hypersetup{
   citecolor=black,
   linkcolor=black,
   filecolor=black,
   urlcolor=blue}
 \else\if@links@bright
  % Use bright colors for the links
  \hypersetup{urlcolor=cyan}
 \else
  % This loads a package that allows extra colors for links.
  \RequirePackage[usenames,dvipsnames]{color}
  % Custom color for references.
  \definecolor{DarkGreen}{rgb}{0,0.6,0}
  % Use more muted colors.
  \hypersetup{
   citecolor=DarkGreen,
   linkcolor=Maroon,
   urlcolor=blue,
   filecolor=black}
 \fi\fi
\fi


%% ----- SET THE FONT SIZES ----------------------------------------------------
% Change the base font size for journal papers.
\if@q@compress
 % 'compress' sizes and skips
 \compresstitle
 \compresssections
 \smallfonts
\else\if@aiaa@journal
 % For journal, set font size to 9.
 \journalfont
\fi\fi



%% ----- BEGIN DOCUMENT AREA ---------------------------------------------------
% List of actions to take after the document has begun.
\AtBeginDocument{ %
 % Commands for each format.
 \if@aiaa@submit %
  % 1.5-spaced lines
 % \onehalfspacing %
 \doublespacing %
 \else\if@aiaa@journal
  % Sloppy hyphenation
  \sloppy
 \fi\fi %
 % Set the PDF title.
 % We have to wait until here to make sure the user has set the title.
 \if@q@links
  \hypersetup{pdftitle=\inserttitle}
 \fi
}
