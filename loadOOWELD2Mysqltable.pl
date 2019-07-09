#!/usr/bin/perl -w

#Load the oo.weld table into mysql table

use strict;
my $description ="Load the oo.weld table into mysql table.\n";

if(scalar(@ARGV) != 3) { die "$description use:\nloadOOWELD2Mysqltable.pl filename dbname table_name\n";}

my $command;
my $filename = $ARGV[0];
my $dbname = $ARGV[1];
my $tablename = $ARGV[2];


$command = "mysql -ulinchen2 -plinchen2 $dbname -e\'
drop table if exists $tablename;
create table $tablename(
QNAME varchar(20),
QB int unsigned not null,
QE int unsigned not null,
QLEN int unsigned not null,
SNAME varchar(20),
SB int unsigned not null,
SE int unsigned not null,
SLEN int unsigned not null,
FRACBPMATCH float,
BPALIGN int unsigned,
SIZEALIGN int unsigned,
SCORE int unsigned,
QDEFN varchar(50),
SDEFN varchar(50),
ERRORB varchar(20),
ERRORE varchar(20),
ROW varchar(80),
FILE varchar(80),
BEGIN int unsigned,
END int unsigned,
indel_N int unsigned,
indel_S int unsigned,
base_S int unsigned,
base_Match int unsigned,
base_Mis int unsigned,
transversions int unsigned,
transitions int unsigned,
per_sim float,
SE_sim float,
per_sim_indel float,
SE_sim_indel float,
K_jc float,
SE_jc float,
k_kimura float,
SE_kimura float,
K_tn float,
SE_tn float,
d_lake float,
AA  int unsigned,
AT int unsigned,
AC int unsigned,
AG int unsigned,
AN int unsigned,
TA int unsigned,
TT int unsigned,
TC int unsigned,
TG int unsigned,
TN int unsigned,
CA int unsigned,
CT int unsigned,
CC int unsigned,
CG int unsigned,
CN int unsigned,
GA int unsigned,
GT int unsigned,
GC int unsigned,
GG int unsigned,
GN int unsigned,
NA int unsigned,
NT int unsigned,
NC int unsigned,
NG int unsigned,
NN int unsigned,
index(QNAME, QB, QE),
index(SNAME, SB, SE),
index(per_sim),
index(base_S));
load data local infile \"$filename\" into table $tablename;
delete from $tablename where qlen=0 or slen = 0;
alter table $tablename add kb int unsigned;
alter table $tablename add k float;
alter table $tablename add type varchar(10);
alter table $tablename add index(kb);
alter table $tablename add index(k);
update $tablename set kb=if(base_s>=50000, 50, if(base_s>=10000, floor( base_s/10000)*10, floor(base_s/1000))), type = if(substring_index(qname, \"_\", 1) = substring_index(sname, \"_\", 1), \"intra\", \"inter\"), k = floor(k_kimura * 10) / 10;
\n\'";

print "$command\n";

system($command);

