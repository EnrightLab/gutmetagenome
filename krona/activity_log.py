#!/usr/bin/env python3 
import pandas as pd
import numpy as np
import dateutil.parser
from datetime import date
from dateutil.parser import *
import datetime
import random
import seaborn as sns
import textwrap
from matplotlib import pyplot as plt
import matplotlib.dates as mdates

fastq=""
with open('/Users/aje/bin/gutmetagenome/krona/fastq.txt', 'r') as file:
    for line in file:
        fastq = fastq + line.strip()

#print(fastq)
        
now = datetime.datetime.now()
current_date_time = now.strftime("%H:%M:%S, %d/%m/%Y")

months=["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]

f = open("/Users/aje/bin/gutmetagenome/krona/cgs.html","w")

def html_h(file,level,text):
    file.write("<H" + str(level) + ">" + text + "</H" + str(level) + ">\n")

def gigabase(total):
    return(total/1000000000)

def megabase(total):
    return(total/1000000)

def terabase(total):
    return(total/1000000000000)

def html_tstart():
    return("<TABLE class=\"styled-table\">")
    
def html_tend():
    return("</TABLE>")
    
def html_trow(name,vals,rr):
    string=""
    string += "<TR>"
    if (name != ""):
        string += str("<TD>" + name + "</TD>")
    for val in vals:
        if (type(val) == float):
            if (rr != 0):
                val=round(val,rr)
            else:
                val=round(val)
        string += str("<TD>" + str(val) + "</TD>",)
    string+= "</TR>"
    return(string)
    
html_blob = textwrap.dedent("""\
<HTML>
<meta http-equiv="refresh" content="60">
<body style="background-color:white;">
<style>
.styled-table {
    width: 30%;
    border-collapse: collapse;
    margin: 25px 0;
    font-size: 1.4em;
    font-family: sans-serif;
    min-width: 400px;
    box-shadow: 0 0 20px rgba(0, 0, 0, 0.15);
}
.styled-table thead tr {
    background-color: #009879;
    color: #ffffff;
    text-align: left;
}
.styled-table th,
.styled-table td {
    padding: 12px 15px;
}
.styled-table tbody tr {
    border-bottom: 1px solid #dddddd;
}

.styled-table tbody tr:nth-of-type(even) {
    background-color: #f3f3f3;
}

.styled-table tbody tr:last-of-type {
    border-bottom: 2px solid #009879;
}
.styled-table tbody tr.active-row {
    font-weight: bold;
    color: #009879;
}

marquee {
  font-size: 80px;
  scrollamount: 30;
  text-align: left;
  font-family: monospace;
}



body {
  font-family: -apple-system, system-ui, Helvetica, Arial, sans-serif;
  margin: 20px;
  display: flex;
  flex-direction: column;
}
</style>
""")

f.write(html_blob)
html_h(f,1,"Genomics Facility Status: " + str(current_date_time))

f.write("<marquee behavior=\"scroll\" scrollamount=\"20\" direction=\"left\" bgcolor=\"lightgrey\">" + fastq + "</marquee>\n")

df = pd.DataFrame(columns=('Platform', 'RunID', 'Date', 'Name', 'Samples', 'Info', 'Length', 'Reads1', 'Reads2', 'Total_bp'))
col_names=('Platform', 'RunID', 'Date', 'Name', 'Samples', 'Info', 'Length', 'Reads1', 'Reads2', 'Total_bp','Age','Recent','Mon','Week','Year','Orig','NiceWeek')


df = pd.read_table('/Volumes/CGS-FS3/SERVICES/Bioinformatic/Anton_Research/sequencing_logs/activity_log.txt',names=col_names) 

monthly_samples=0
monthly_reads=0
monthly_len=0
monthly_bp=0

all_samples=0
all_reads=0
all_len=0
all_bp=0

df["Orig"]=df.index

for index, row in df.iterrows():
    thisday=dateutil.parser.parse(row['Date'],dayfirst = True)
    delta=(now-thisday).days
    month_yr= str(thisday.month)

    month_int= int(thisday.month)
    month_name = months[month_int-1]
    week_yr= str(thisday.isocalendar()[1])
    yr=(str(thisday.year))

    df = df.astype({'NiceWeek': str})

    df.iloc[index,df.columns.get_loc('NiceWeek')]=str(yr + " " + month_name)
    df.iloc[index,df.columns.get_loc('Mon')]=int(month_yr)
    df.iloc[index,df.columns.get_loc('Week')]=int(week_yr)
    df.iloc[index,df.columns.get_loc('Year')]=int(yr)
    df.iloc[index,df.columns.get_loc('Age')]=delta
    
    if(delta <= 31):
        #print(row['Platform'],row['RunID'],delta,row['Info'])
        monthly_samples+=row['Samples']
        monthly_reads+=row['Reads2']
        monthly_len+=row['Length']
        monthly_bp+=row['Total_bp']
        df.iloc[index,df.columns.get_loc('Recent')]=1
    else:
        df.iloc[index,df.columns.get_loc('Recent')]=0

weekly=df[df["Age"] <= 7]
monthly=df[df["Age"] <= 31]
yearly=df[df["Age"] <= 365]

df["Experiments"]=1

summary=np.empty([3,12])

platforms=["MISEQ","NEXSEQ","PROMETHION"]
index=0
for platform in (platforms):
    dfw=weekly[weekly["Platform"] == platform]
    dfm=monthly[monthly["Platform"] == platform]
    dfy=yearly[yearly["Platform"] == platform]
    dfa=df[df["Platform"] == platform]
    
    summary[index,0]=int(dfw["Samples"].sum())
    summary[index,1]=int(dfm["Samples"].sum())
    summary[index,2]=int(dfy["Samples"].sum())
    summary[index,3]=int(dfa["Samples"].sum())

    summary[index,4]=gigabase(int(dfw["Total_bp"].sum()))
    summary[index,5]=gigabase(int(dfm["Total_bp"].sum()))
    summary[index,6]=gigabase(int(dfy["Total_bp"].sum()))
    summary[index,7]=gigabase(int(dfa["Total_bp"].sum()))

    summary[index,8]=int(dfw["Samples"].count())
    summary[index,9]=int(dfm["Samples"].count())
    summary[index,10]=int(dfy["Samples"].count())
    summary[index,11]=int(dfa["Samples"].count())
    
    index=index+1


summary = pd.DataFrame(summary,dtype=str,index=pd.Index(['Illumina MiSeq', 'Illumina NextSeq', 'Nanopore PromethION'], name='Rows'),columns=pd.Index(['Weekly Samples', 'Monthly Samples', 'Yearly Samples', 'All Samples','Weekly Bp', 'Monthly Bp', 'Yearly Bp', 'All Bp','Weekly Ex', 'Monthly Ex', 'Yearly Ex', 'All Ex'], name='Cols'))
summary = summary.astype({"Weekly Bp": float, "Monthly Bp": float, "Yearly Bp": float, "All Bp": float})


exp = summary.iloc[:, 8:12]
samples = summary.iloc[:, 0:4].copy()
bases = summary.iloc[:, 4:8]

i=0
for index,row in samples.iterrows():
    for j in range(len(row)):
        samples.iloc[int(i),int(j)]=str(round(float(samples.iloc[int(i),int(j)]))) + " (" + str(round(float(exp.iloc[int(i),int(j)]))) + ")"
    i=i+1

f.write("<TABLE BORDER=0><TR><TD>")
f.write(html_tstart())
html_h(f,3,"GigaBases Sequenced")
header=["<b>Platform</b>","<b>Week</b>","<b>Month</b>","<b>Year</b>","<b>Total</b>"]
f.write(html_trow("",header,2))

for index, row in bases.iterrows():
    f.write(html_trow(row.name,row,2))
f.write(html_tend())
f.write("</TD><TD><IMG SRC=\"fig1.png\"><IMG SRC=\"fig2.png\">")
f.write("</TD></TR>")
f.write("<TR><TD>")

html_h(f,3,"No. Samples (No. Experiments)")
f.write(html_tstart())
header=["<b>Platform</b>","<b>Week</b>","<b>Month</b>","<b>Year</b>","<b>Total</b>"]
f.write(html_trow("",header,2))
for index, row in samples.iterrows():
    f.write(html_trow(row.name,row,0))
f.write(html_tend())
f.write("</TD><TD>")
#html_h(f,4,"Experiments")
#f.write(html_tstart())
#header=["Platform","This Week","This Month","This Year","Total"]
#f.write(html_trow("",header,2))
#for index, row in exp.iterrows():
#    f.write(html_trow(row.name,row,0))
#f.write(html_tend())
f.write("<IMG SRC=\"fig4.png\">")
f.write("</TD></TR>")
f.write("</TABLE>")



df['Total_bp'] = df['Total_bp'].div(1000000000).round(2)
df_new=df[df["Age"] <= (12*7)]

#display(df_new)
df2=df_new.groupby(["Year","Mon","Week","NiceWeek","Platform"])["Total_bp"].sum().unstack(fill_value=0).stack().reset_index()
df2.columns=["Year","Mon","Week","NiceWeek","Platform","Sum"]
df2=df2.sort_values(["Year","Week"],ascending=[True,True])

pivot=pd.pivot_table(df_new,index=['Year','Mon','Week','NiceWeek'], columns=['Platform'],aggfunc="sum", fill_value=0).reset_index()


width = 0.8
plt.figure(figsize=(6.6, 4))
plt.bar(pivot.index,pivot["Total_bp"]["PROMETHION"]+pivot["Total_bp"]["NEXSEQ"]+pivot["Total_bp"]["MISEQ"],label="Illumina MiSeq",color='orange')
plt.bar(pivot.index,pivot["Total_bp"]["PROMETHION"]+pivot["Total_bp"]["NEXSEQ"],label="Illumina NextSeq",color='red')
plt.bar(pivot.index,pivot["Total_bp"]["PROMETHION"],label="Nanopore",color='darkblue')
plt.xticks(rotation = 90)
plt.xticks(pivot.index,pivot["NiceWeek"].values)
plt.title('Nucleotides Sequenced (Past 3 months)')
plt.xlabel("Week")
plt.ylabel("Total (Gigabases)")
plt.legend()
plt.subplots_adjust(top=0.925, bottom=0.20, left=0.09, right=0.99, hspace=0.01, wspace=0.01)
plt.savefig('fig1.png')
#plt.show()
plt.close

width = 0.8
plt.figure(figsize=(6.6, 4))
plt.bar(pivot.index,pivot["Samples"]["PROMETHION"]+pivot["Samples"]["NEXSEQ"]+pivot["Samples"]["MISEQ"],label="MiSeq",color='orange')
plt.bar(pivot.index,pivot["Samples"]["PROMETHION"]+pivot["Samples"]["NEXSEQ"],label="NextSeq",color='red')
plt.bar(pivot.index,pivot["Samples"]["PROMETHION"],label="Nanopore",color='darkblue')
plt.xticks(rotation = 90)
plt.xticks(pivot.index,pivot["NiceWeek"].values)
plt.title('No of Samples (Past 3 months)')
plt.xlabel("Week")
plt.ylabel("Total")
plt.legend()
plt.subplots_adjust(top=0.925, bottom=0.20, left=0.09, right=0.99, hspace=0.01, wspace=0.01)
plt.savefig('fig2.png')
#plt.show()


pivot=pd.pivot_table(df_new,index=['Year','Mon','Week','NiceWeek'], columns=['Platform'],aggfunc=len, fill_value=0).reset_index()

#width = 0.8
#plt.figure(figsize=(10, 4))
#plt.bar(pivot.index,pivot["Total_bp"]["PROMETHION"]+pivot["Total_bp"]["NEXSEQ"]+pivot["Total_bp"]["MISEQ"],label="MiSeq",color='orange')
#plt.bar(pivot.index,pivot["Total_bp"]["PROMETHION"]+pivot["Total_bp"]["NEXSEQ"],label="NextSeq",color='red')
#plt.bar(pivot.index,pivot["Total_bp"]["PROMETHION"],label="Nanopore",color='darkblue')
#plt.xticks(rotation = 90)
#plt.xticks(pivot.index,pivot["NiceWeek"].values)
#plt.title('No of Experiments')
#plt.legend()
#plt.xlabel("Week")
#plt.ylabel("Total")
#plt.savefig('fig3.png')
#plt.show()




df_new=df
pivot=pd.pivot_table(df_new,index=['Year','Mon','Week','NiceWeek'], columns=['Platform'],aggfunc=len, fill_value=0).reset_index()

pp=pd.DataFrame(pivot["Total_bp"]["PROMETHION"] + pivot["Total_bp"]["MISEQ"] + pivot["Total_bp"]["NEXSEQ"],columns=["Experiments"])
pivot=pd.pivot_table(df_new,index=['Year','Mon','Week','NiceWeek'], columns=['Platform'],aggfunc="sum", fill_value=0).reset_index()

pp["Total_bp"]=pivot["Total_bp"]["PROMETHION"] + pivot["Total_bp"]["MISEQ"] + pivot["Total_bp"]["NEXSEQ"]
pp["Samples"]=pivot["Samples"]["PROMETHION"] + pivot["Samples"]["MISEQ"] + pivot["Samples"]["NEXSEQ"]
pp["Week"]=pivot["NiceWeek"]

pp["Experiments Avg"]=pp["Experiments"].rolling(4,center=True).mean()
pp["Experiments Avg"]=pp["Experiments Avg"]
pp["Bp Avg"]=pp["Total_bp"].rolling(4,center=True).mean()
pp["Samples Avg"]=pp["Samples"].rolling(4,center=True).mean()

#display(pp)

plt.figure(figsize=(13.7, 5.0))
sns.set_theme(rc={'figure.figsize':(14.7,5.0)})
plt=sns.lineplot(x="Week",y="Samples Avg",data=pp,color="orange",label="No. Samples")
plt=sns.lineplot(x="Week",y="Bp Avg",data=pp,color="red",label="Total Gigabases")
plt.set(ylabel='Total Samples\nTotal Bases (Gigabases)')

sns.move_legend(plt, "upper left")
plt.tick_params(axis='x', labelrotation = 90)
plt.grid(None) 
ax2 = plt.twinx()
ax2.set(ylim=(0, 40))
plt=sns.lineplot(x="Week",y="Experiments Avg",data=pp,color="darkblue",label="No. Experiments",errorbar=None)
plt.set_title("4 Week Averages")
plt.grid(None) 

fig = plt.get_figure()

fig.subplots_adjust(top=0.925, bottom=0.20, left=0.060, right=0.90, hspace=0.01, wspace=0.01)
fig.savefig("fig4.png")
#f.write("<IMG SRC=\"fig4.png\">")

f.write("</HTML>")
f.flush
f.close()
