#!/usr/bin/perl -w
# perl -w in case I mess the line endings up
use strict;
use warnings;

# SETTINGS ############################################
my $indent = "  "; # <--- SET INDENTATION SEQUENCE HERE
#######################################################



#############################################################################
##
##    Outline Beamer Class Presentation Maker
##    Copyright (C 2008 Jan Schejbal
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 2 of the License, or
##    (at your option) version 3 of the License.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>
##    or write to the Free Software Foundation, Inc.,
##    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
##
###############################################################################








my $currlevel = 0;

## THESE ARRAYS/BLOCKS ARE IN THE PUBLIC DOMAIN ###############################
# As parts of these arrays get copied into the output .tex, they are PD so
# the users of this tool keep full copyright on the output without any doubt
my @itemstart = (
  '\\section{',
  '  \\subsection{',
  '    \\begin{frame}[fragile] \\frametitle{', # fragile to allow easy \verbose
  '        \\item ',
  '          \\item ',
  '            \\item '
);

my @itemclose = (
  "}\n",
  "}\n",
  "}\n",
  "\n",
  "\n",
  "\n"
);

my @itemend = (
  "\n\n\n",
  "\n\n\n",
  "\\end{frame}\n\n", # indentation NOT allowed because of [fragile] option
  "",
  "",
  ""
);

my @prelevel = (
  "",
  "",
  "",
  "      \\begin{itemize}\n",
  "        \\begin{itemize}\n",
  "          \\begin{itemize}\n"
);

my @postlevel = (
  "",
  "",
  "",
  "      \\end{itemize}\n",
  "        \\end{itemize}\n",
  "          \\end{itemize}\n"
);


print <<END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    This LaTeX source was automatically generated by outlinebeamer.pl     %
%                                                                          %
%   Remember that most changes should be made to the outline source file   %
% changes to this file will be overwritten when the file gets regenerated  %
%                                                                          %
%                  outlinebeamer - Outline Beamer Class Presentation Maker %
%                                    http://outlinebeamer.sourceforge.net/ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END
;
## END OF PUBLIC DOMAIN BLOCK ################################################

my $error = 0;
my %config;

$config{quotes} = 0;
$config{passleft}  = 0;
$config{subtractmaxindent} = 0;
$config{blocks} = 0;
my $openblocks = 0;
my @blockstack = ();

my $passingtrough = 0;
my $endtag = "";


while (<>) {
  if ( $passingtrough ) { # handle passtrough mode
    if ( m/^$endtag>>/ ) {
      $passingtrough = 0;
      next;
    } else {
      print;
      next;
    }
  }
  if ( !s/^@@//) {  # if line starts with @@, remove that and skip tag interpretation
    if ( m/^@<<(\w*)/ ) { # if line starts with @<<, enable passtrough mode, use the optional word after << as end tag
      $passingtrough = 1;
      $endtag = $1;
      next;
    }
    if ( m/^#/ ) { next; } # if line starts with #, ignore it
    if ( s/^@!//) { print; next;} # if line starts with @!, pass it unchanged (after removing @!)
    if ( m/^@\+(\w*)/) { # if line starts with @+xxxx, enable feature xxxx, ignore line
      if ( ! defined $config{$1} ) {
        print STDERR "Warning: unsuppored option $1 used at line $.\n";
      }
      $config{$1} = 1;
      next;
    }
    if ( s/^@\-(.*)\s*//) { # if line starts with @-xxxx, disable feature xxxx, ignore line
      if ( $config{$1} != 1 ) {
        print STDERR "Warning: attempted to disable non-enabled option $1 at line $.\n";
      }
      $config{$1} = 0;
      next;
    }
    if ( m/^\@end/ ) { # if the line starts with @end, end outline
      while ( $currlevel > 0 ) {
        print $itemend[$currlevel];
        print $postlevel[$currlevel];
        $currlevel--;
      }
      print $itemend[0];
      next;
    }
  }  
  
  if ( m/^\s*$/ ) { next; }  # skip empty lines
  s/\s*$//; # clean white space at end of line - for passing on lines insert \n!
  if ( $config{quotes} ) {
    s/(\W)"(\w)/$1\\glqq{}$2/g;
    s/(\w)"(\W)/$1\\grqq{}$2/g;
  }
  my $thislevel = 0;
  my $oldlevel = $currlevel;
  
  while (s/^$indent//) { $thislevel++; }
  if ( $config{passleft} ) { # if passleft ("zero indent = plain") option enabled, correct indentation count
    $thislevel--;
    if ( $thislevel == -1 ) {   # pass zero indent lines unchanged
      print; print "\n"; 
      next;
    }
  }
  if ( $thislevel > 5 ) {  # pass max indent lines
    print $indent x ($config{subtractmaxindent} ? $thislevel - 6 : $thislevel); # keep indentation
    print; print "\n"; 
    next;
  }
  
  if ( $config{blocks} && $thislevel == 3 ) {
    unless ( s/^@@\[/[/ || s/^@@\]/]/ ) { # if blocks enabled, allow @@ after indentation to disable block feature
      if ( m/^\[([a-zA-Z]*)\[(.*)\s*$/ ) {
        my $envname = ( $1 eq "" ) ? "block" : $1; # environment defaults to block
        my $param = ( $2 ne "" ) ? "{$2}" : ""; # create parameter brackets only if parameter given
        print "      \\begin{$envname}$param\n";
        $blockstack[$openblocks++] = $envname;
        next;
      }
      if ( m/^\]\]/ ) {
        if ( $openblocks > 0 ) {
          print "      \\end{".$blockstack[--$openblocks]."}\n";
        } else {
          print STDERR "ERROR: Attempted to close block while no block open\n";
          $error = 1;
        }
        next;
      }
    }
  }
  
  while ( $currlevel < $thislevel ) {
    $currlevel++;
    print $prelevel[$currlevel];
  }
  while ( $currlevel > $thislevel ) {
    print $itemend[$currlevel];
    print $postlevel[$currlevel];
    $currlevel--;
  }
  if ( $thislevel <= $oldlevel ) {
    print $itemend[$thislevel];
  }
  
  print $itemstart[$thislevel].$_.$itemclose[$thislevel];
  
  if ( $thislevel > 2 && $oldlevel < 2 ) {
    print STDERR "ERROR: Starting item without frame at line $., tex will probably not compile\n";
    $error = 1;
  }
  
}

if ( $currlevel > 0 ) {
  print STDERR "ERROR: File not ending on level 0 - forgot \@end or garbage after \@end?\n";
  $error = 1;
}

if ( $openblocks > 0 ) {
  print STDERR "ERROR: Not all blocks were closed - forgot ']]'?\n";
  $error = 1;
}

exit($error);

