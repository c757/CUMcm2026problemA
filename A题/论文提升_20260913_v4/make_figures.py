from pathlib import Path
import json
import logging
import os
import sys
from solve import ROOT, properties
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'.work/mplconfig'))
os.environ.setdefault('SOURCE_DATE_EPOCH','1788998400')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter
from utils.setup_style import setup_style, CJK_FONT_PRIORITY
from utils.plot_style import export_figure as audited_export
from utils.export_figure import export_figure as pdf_export

OUT=ROOT/'figures'
DATA=ROOT/'results'
BLUE='#0072B2'; ORANGE='#D55E00'; GREEN='#009E73'; GRAY='#666666'
COLORS=[BLUE,ORANGE,GREEN]
LINES=['-','--','-.']
WIDTH=6.3

class GlyphFailureHandler(logging.Handler):
    def emit(self,record):
        message=record.getMessage()
        if 'does not have a glyph' in message or 'missing from font' in message:
            raise RuntimeError(message)

def configure():
    logging.getLogger('matplotlib').addHandler(GlyphFailureHandler())
    sources=[ROOT/'assets/wqy-zenhei.ttc',
             ROOT/'.work/tex_runtime/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc',
             Path('/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc')]
    source=next((p for p in sources if p.exists()),None)
    if source is None:
        raise RuntimeError('Install fonts-wqy-zenhei or supply assets/wqy-zenhei.ttc')
    font_manager.fontManager.addfont(str(source))
    CJK_FONT_PRIORITY.insert(0,'WenQuanYi Zen Hei')
    setup_style(journal='general',lang='zh',use_sciplots=False)
    plt.rcParams.update({'font.size':8.5,'axes.labelsize':9,'axes.titlesize':9,
                         'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,
                         'axes.titlelocation':'left','axes.spines.top':False,
                         'axes.spines.right':False,'lines.linewidth':1.4,
                         'figure.facecolor':'white','axes.facecolor':'white',
                         'legend.frameon':False,'savefig.dpi':300,'svg.hashsalt':'20260910'})
    OUT.mkdir(exist_ok=True)

def panels(count=1,height=2.9):
    return plt.subplots(1,count,figsize=(WIDTH,height),layout='constrained',squeeze=False)[0:2]

def save(fig,stem):
    for ax in fig.axes:
        ax.minorticks_off()
        if ax.get_yscale()=='log' and not hasattr(ax,'_colorbar'):
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:.0e}'))
        if hasattr(ax,'_colorbar') and ax._colorbar.solids is not None:
            ax._colorbar.solids.set_rasterized(False)
    audited_export(fig,OUT/stem,dpi=300,strict_layout=True,strict_design=True)
    pdf_export(fig,str(OUT/stem),formats=['pdf'],dpi=300,
               size_inches=tuple(fig.get_size_inches()),grayscale_preview=False,tight=False)
    plt.close(fig)

def series(ax,x,ys,labels):
    for y,label,color,ls in zip(ys,labels,COLORS,LINES):
        ax.plot(x,y,label=label,color=color,ls=ls)
    ax.legend()
    ax.margins(x=0.01)

def make_figures():
    configure()
    env=pd.read_csv(DATA/'input_environment.csv')
    rad=pd.read_csv(DATA/'input_radius.csv')
    summary=json.loads((DATA/'summary.json').read_text())
    q1=np.load(DATA/'q1_fields.npz');q3=np.load(DATA/'q3_fields.npz');q4=np.load(DATA/'q4_fields.npz')

    fig,axes=panels(2);a,b=axes[0]
    a.plot(env.time_s/3600,env.temperature_C,color=BLUE)
    b.plot(env.time_s/3600,env.moisture_kg_kg,color=ORANGE)
    for ax in (a,b):ax.set_xlabel('时间 / h');ax.set_xlim(0,4)
    a.set_ylabel('环境温度 / ℃');b.set_ylabel('有效平衡含水率 / (kg/kg)')
    a.set_title('(a) 温度边界');b.set_title('(b) 水分边界')
    save(fig,'raw_q1_environment')

    tail=env[env.time_s>=10800]
    fig,axes=panels(2);a,b=axes[0]
    a.hist(tail.temperature_C,bins=8,color=BLUE,alpha=.75,edgecolor='white')
    b.hist(tail.moisture_kg_kg,bins=8,color=ORANGE,alpha=.75,edgecolor='white')
    a.axvline(50,color='black',ls='--',lw=1);b.axvline(.05,color='black',ls='--',lw=1)
    a.set_xlabel('环境温度 / ℃');b.set_xlabel('有效平衡含水率 / (kg/kg)')
    for ax in (a,b):ax.set_ylabel('频数');ax.set_title('末小时：61个记录')
    b.ticklabel_format(axis='x',style='plain',useOffset=False)
    b.set_xticks([.0498,.05,.0502])
    save(fig,'raw_q2_tail')

    fig,axes=panels(height=2.7);a=axes[0,0];c=np.linspace(.05,2.55,501)
    for t,ls,col in ((28,'--',BLUE),(50,'-',ORANGE)):
        a.semilogy(c,properties(3,np.full_like(c,t),c)[2],color=col,ls=ls,label=f'{t}℃')
    a.set_xlabel('干基含水率 C / (kg/kg)');a.set_ylabel('扩散系数 D / (m²/s)')
    a.legend();a.set_title('题设物性函数（非实测）')
    save(fig,'raw_q3_diffusivity')

    fig,axes=panels(height=2.7);a=axes[0,0]
    a.plot(rad.time_s/3600,rad.radius_cm,color=BLUE,lw=1)
    a.scatter(rad.time_s/3600,rad.radius_cm,s=7,color=BLUE)
    a.set_xlabel('时间 / h');a.set_ylabel('药材半径 / cm');a.set_xlim(0,72)
    a.set_title('题给半径：145个记录')
    save(fig,'raw_q4_radius')

    fig,axes=panels(height=2.7);a=axes[0,0]
    times=q1['time_s']
    series(a,times/60,[q1['Tcenter'],q1['Tsurface'],np.interp(times,env.time_s,env.temperature_C)],
           ['中心','表面','环境'])
    a.set_xlabel('时间 / min');a.set_ylabel('温度 / ℃')
    save(fig,'process_q1_heating')

    fig,axes=panels(2);a,b=axes[0]
    for t,color,ls in zip([100,600,1800],COLORS,LINES):
        i=np.flatnonzero(times==t)[0]
        a.plot(q1['x']*2,q1['temperature'][i],color=color,ls=ls,label=f'{t} s')
        b.plot(q1['x']*2,q1['moisture'][i],color=color,ls=ls,label=f'{t} s')
    a.set_ylabel('温度 / ℃');b.set_ylabel('干基含水率 / (kg/kg)')
    for ax in (a,b):ax.set_xlabel('半径位置 / cm');ax.legend()
    a.set_title('(a) 温度');b.set_title('(b) 水分')
    save(fig,'result_q1_profiles')

    take=q3['time_s']<=10800;tt=q3['time_s'][take]/3600;xx=q3['x']*2
    temps=q3['temperature'][take];moist=q3['moisture'][take]
    diff=properties(3,temps,moist)[2]
    fig,axes=panels(height=2.9);a=axes[0,0]
    mesh=a.pcolormesh(tt,xx,diff.T,norm=LogNorm(),cmap='viridis',shading='auto')
    cb=fig.colorbar(mesh,ax=a,ticks=np.geomspace(diff.min(),diff.max(),4),
                    format=FuncFormatter(lambda x,pos:f'{x:.1e}'))
    cb.minorticks_off();cb.set_label('D / (m²/s)，对数色标')
    a.set_xlabel('时间 / h');a.set_ylabel('半径位置 / cm')
    save(fig,'process_q2_diffusion')

    fig,axes=panels(2,height=3.1);a,b=axes[0]
    for ax,field,label in [(a,temps,'温度 / ℃'),(b,moist,'C / (kg/kg)')]:
        mesh=ax.pcolormesh(tt,xx,field.T,cmap='viridis',shading='auto')
        cb=fig.colorbar(mesh,ax=ax,pad=.02);cb.set_label(label)
        ax.set_xlabel('时间 / h');ax.set_ylabel('半径位置 / cm')
    a.set_title('(a) 温度场');b.set_title('(b) 含水率场')
    save(fig,'result_q2_fields')

    fig,axes=panels(height=2.9);a=axes[0,0]
    series(a,q3['time_s']/3600,[q3['Cmax'],q3['Cmean'],q3['Csurface']],['全域最大','径向体积权均值','表面'])
    a.axhline(.15,color='black',ls=':',lw=1.1)
    a.set_xlabel('时间 / h');a.set_ylabel('干基含水率 / (kg/kg)')
    save(fig,'process_q3_threshold')

    fig,axes=panels(height=2.7);a=axes[0,0]
    a.plot(xx,q3['moisture'][-1],color=BLUE,label='严格终点的含水率')
    a.axhline(.15,color='black',ls='--',label='阈值 0.15')
    a.set_xlabel('半径位置 / cm');a.set_ylabel('干基含水率 / (kg/kg)')
    a.legend();a.set_ylim(0,.17)
    save(fig,'result_q3_terminal')

    fig,axes=panels(height=3.1);a=axes[0,0]
    ts=q4['time_s']/3600;rs=q4['radius_cm'];rr=np.linspace(0,2,161)
    physical=np.full((len(ts),len(rr)),np.nan)
    for i,radius in enumerate(rs):
        inside=rr<=radius
        physical[i,inside]=np.interp(rr[inside]/radius,q4['x'],q4['moisture'][i])
    mesh=a.contourf(ts,rr,np.ma.masked_invalid(physical.T),cmap='viridis',
                    levels=np.linspace(0,2.55,22),corner_mask=False)
    fig.colorbar(mesh,ax=a).set_label('C / (kg/kg)')
    a.plot(ts,rs,color='black',lw=1,label='实际表面')
    a.set_xlabel('时间 / h');a.set_ylabel('物理半径位置 / cm');a.legend(loc='upper right')
    save(fig,'process_q4_moving')

    fig,axes=panels(height=2.7);a=axes[0,0]
    labels=['问题3：固定半径','问题4物性：固定半径','问题4：给定收缩']
    hours=[summary['runs']['3']['crossing_h'],summary['counterfactual']['crossing_h'],
           summary['runs']['4']['crossing_h']]
    for i,(value,color) in enumerate(zip(hours,COLORS)):
        a.scatter(value,i,color=color,s=32,zorder=3)
        a.annotate(f'{value:.3f} h',(value,i),xytext=(6,5),textcoords='offset points',fontsize=8)
    a.set_yticks(range(3),labels);a.set_xlabel('全域达到阈值的事件时间 / h')
    a.set_xlim(0,max(hours)*1.18);a.set_ylim(2.5,-.5);a.grid(axis='x',alpha=.2)
    save(fig,'result_q4_counterfactual')

    grid=pd.read_csv(DATA/'grid_comparison.csv')
    fig,axes=panels(2);a,b=axes[0]
    for q,color,ls in zip([1,2,3,4],[BLUE,ORANGE,GREEN,GRAY],['-','--','-.',':']):
        frame=grid[grid.question==q]
        a.loglog(frame.fine_n,frame.moisture_max_difference,color=color,ls=ls,marker='o',label=f'问题{q}')
        b.loglog(frame.fine_n,frame.temperature_max_difference,color=color,ls=ls,marker='o',label=f'问题{q}')
    a.axhline(5e-5,color='black',lw=.8);b.axhline(.005,color='black',lw=.8)
    for ax in (a,b):
        ax.set_xlabel('细网格区间数 N');ax.legend(ncol=2)
        ticks=sorted(grid.fine_n.unique());ax.set_xticks(ticks,[str(t) for t in ticks])
    a.set_ylabel('含水率最大差 / (kg/kg)');b.set_ylabel('温度最大差 / ℃')
    save(fig,'process_q3_convergence')

    sens=pd.read_csv(DATA/'sensitivity.csv')
    names=['last','mean','T_minus_0.5','T_plus_0.5','Ce_minus_0.001','Ce_plus_0.001','hm_minus_10pct','hm_plus_10pct']
    labels=['末值延拓','末小时均值延拓','尾段温度 -0.5℃','尾段温度 +0.5℃','尾段平衡值 -0.001','尾段平衡值 +0.001','传质系数 -10%','传质系数 +10%']
    fig,axes=panels(height=3.2);a=axes[0,0]
    for q,col,offset,marker in [(3,BLUE,-.12,'o'),(4,ORANGE,.12,'s')]:
        frame=sens[sens.question==q].set_index('scenario')
        vals=frame.loc[names,'change_pct'].values
        a.scatter(vals,np.arange(8)+offset,color=col,marker=marker,s=22,label=f'问题{q}')
    a.set_yticks(range(8),labels);a.invert_yaxis();a.axvline(0,color='black',lw=.7)
    a.set_xlabel('事件时间相对名义情景的变化 / %');a.legend(loc='lower right')
    save(fig,'result_q4_sensitivity')
    print(json.dumps({'figures':len(list(OUT.glob('*.svg'))),'format':['svg','png','pdf'],'dpi':300}))

if __name__=='__main__':make_figures()





