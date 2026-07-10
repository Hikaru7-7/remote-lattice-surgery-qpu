import json
import sys
MODE = sys.argv[1]
import os
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, f"chip_frames_{MODE}.json")))
G = d["geom"]; FR = d["frames"]; IONS = d["ions"]
CAPMAP = {}
for x in G["MX"]: CAPMAP[round(x,1)] = 2
for x in G["WELL"]: CAPMAP[round(x,1)] = 4
CAPMAP[round(G["SWX"],1)] = 3
for x in G["PARKX"]: CAPMAP[round(x,1)] = 2
for x in G["SPX"]: CAPMAP[round(x,1)] = 2
CAPMAP[round(G["CAVX"],1)] = 2; CAPMAP[round(G["YBX"],1)] = 1

html = """<!DOCTYPE html><html><head><meta charset="utf-8">
<title>d=3 merge on the chapter-4 chip geometry</title><style>
body{background:#232320;color:#ddd;font:14px -apple-system,'Segoe UI',sans-serif;margin:0;padding:14px}
#hdr{margin:4px 8px 10px}#hdr b{font-size:16px}
#stage{position:relative;background:#232320;border:1px solid #3a3a36;border-radius:10px;overflow-x:auto}
#world{position:relative;width:1330px;height:640px}
.zone{position:absolute;border-radius:8px;border:1px dashed #4c4c46}
.zlbl{position:absolute;color:#8f8f86;font-size:12px}
.well{position:absolute;width:34px;height:34px;border:1.4px solid #56564e;border-radius:9px;transform:translate(-50%,-50%)}
.well.gate{border-color:#4a6f96}.well.swapw{border-color:#7d74c9;border-style:dashed}
.well.spam{border-color:#a2793c}.well.cav{border-color:#3f8f71}.well.park{border-color:#777}
.jcol{position:absolute;width:0;border-left:2px dotted #5c5c54;transform:translateX(-50%)}
.jlink{position:absolute;width:2px;background:#4c4c46;transform:translateX(-50%)}
.jdot{position:absolute;width:7px;height:7px;border-radius:50%;background:#8a8a80;transform:translate(-50%,-50%)}
.wall{position:absolute;width:7px;background:#6e6e66;border-radius:2px;transform:translateX(-50%)}
.ion{position:absolute;width:27px;height:27px;border-radius:8px;display:flex;align-items:center;justify-content:center;
 font-weight:600;font-size:11.5px;transform:translate(-50%,-50%);transition:left .24s ease,top .24s ease;z-index:5;color:#fff}
.data{background:#2e2e2a;border:1.6px solid #9a9a90;color:#e8e8e0}
.X{background:#3f7fd4;border:1.6px solid #6fa5ec}.Z{background:#c96f3b;border:1.6px solid #e59a6c}
.comm{background:#2f8f6b;border:1.6px solid #5cc39a}.spare{background:#3a3a36;border:1.6px solid #6e6e66;color:#aaa}
.yb{width:16px;height:16px;border-radius:50%;background:#d85a30;border:1.4px solid #f08a60;font-size:0}
.hi{box-shadow:0 0 0 3px rgba(240,220,120,.85),0 0 14px rgba(240,220,120,.5);z-index:8}
#cap{min-height:44px;margin:12px 8px 6px;font-size:15px;line-height:1.45}
#badge{display:inline-block;background:#3a3a2f;border:1px solid #6b6b4f;color:#e6d98a;border-radius:20px;
 padding:2px 12px;font-size:12px;margin-left:10px;vertical-align:2px}
#verify{margin:0 8px;font-size:13px}.ok{color:#7ec98f}.bad{color:#e8806a}
#bar{display:flex;gap:10px;align-items:center;margin:10px 8px}
button{background:#32322e;color:#ddd;border:1px solid #55554d;border-radius:8px;padding:6px 14px;font-size:14px;cursor:pointer}
button:hover{background:#3d3d38}input[type=range]{flex:1}
.legend{margin:8px;color:#9a9a90;font-size:12.5px}
.lg{display:inline-block;width:13px;height:13px;border-radius:4px;vertical-align:-2px;margin:0 4px 0 12px}
</style></head><body>
<div id="hdr"><b>__TITLE__</b><br>
<span style="color:#9a9a90">Frames come from the released scheduler's op stream. A strict verifier runs live on every frame:
well capacity, nothing at rest on a junction column, and no reordering along a row without a shared-well crystal rotation.</span></div>
<div id="stage"><div id="world"></div></div>
<div id="cap"></div><div id="verify"></div>
<div id="bar"><button id="b0">&#9198;</button><button id="bp">&#9664;</button><button id="pl">Play</button>
<button id="bn">&#9654;</button><button id="be">&#9197;</button>
<input type="range" id="sl" min="0" value="0"><span id="ctr"></span>
<select id="sp"><option value="900">slow</option><option value="450" selected>normal</option><option value="180">fast</option></select></div>
<div class="legend"><span class="lg" style="background:#2e2e2a;border:1.6px solid #9a9a90"></span>data
<span class="lg" style="background:#3f7fd4"></span>X check <span class="lg" style="background:#c96f3b"></span>Z check
<span class="lg" style="background:#2f8f6b"></span>comm Ba <span class="lg" style="background:#d85a30;border-radius:50%"></span>Yb coolant
<span class="lg" style="border:1.4px dashed #7d74c9;background:none"></span>swap well
<span class="lg" style="border-left:2px dotted #5c5c54;background:none;width:2px"></span>junction column &middot;
memory homes 2d+1 &middot; gate strip: junction, well, junction, well, junction, well, swap &middot; SPAM d sites &middot; wall &middot; cavity</div>
<script>
const G=__GEOM__, FR=__FRAMES__, KIND=__IONS__, CAP=__CAP__;
const W=document.getElementById('world');
function zone(x1,x2,y,cls,lbl){const z=document.createElement('div');z.className='zone';
 z.style.left=(x1)+'px';z.style.width=(x2-x1)+'px';z.style.top=(y-34)+'px';z.style.height='68px';W.appendChild(z);
 if(lbl){const t=document.createElement('div');t.className='zlbl';t.textContent=lbl;t.style.left=x1+'px';t.style.top=(y-52)+'px';W.appendChild(t);}}
const CY=G.CY;
for(const r in CY){const y=CY[r];
 zone(G.MX[0]-26,G.MX[6]+26,y,'', r=='0'?'Memory (2d+1 homes)':'');
 zone(G.JCOL[0]-22,G.WELL[2]+24,y,'', r=='0'?'Gate strip (j w j w j w)':'');
 zone(G.SWX-24,G.SWX+24,y,'', r=='0'?'swap':'');
 zone(G.SPX[0]-26,G.SPX[2]+26,y,'', r=='0'?'SPAM (d sites)':'');
 zone(G.CAVX-26,G.YBX+22,y,'', r=='0'?'Optical I/F':'');
 for(const x of G.MX){const w=document.createElement('div');w.className='well';w.style.left=x+'px';w.style.top=y+'px';W.appendChild(w);}
 for(const x of G.WELL){const w=document.createElement('div');w.className='well gate';w.style.left=x+'px';w.style.top=y+'px';W.appendChild(w);}
 for(const x of G.JCOL){const j=document.createElement('div');j.className='jcol';j.style.left=x+'px';j.style.top=(y-26)+'px';j.style.height='52px';W.appendChild(j);}
 {const w=document.createElement('div');w.className='well swapw';w.style.left=G.SWX+'px';w.style.top=y+'px';W.appendChild(w);}
 if(r=='2'){for(const x of G.PARKX){const w=document.createElement('div');w.className='well park';w.style.left=x+'px';w.style.top=y+'px';W.appendChild(w);}}
 for(const x of G.SPX){const w=document.createElement('div');w.className='well spam';w.style.left=x+'px';w.style.top=y+'px';W.appendChild(w);}
 const wl=document.createElement('div');wl.className='wall';wl.style.left=G.WALLX+'px';wl.style.top=(y-30)+'px';wl.style.height='60px';W.appendChild(wl);
 const cv=document.createElement('div');cv.className='well cav';cv.style.left=G.CAVX+'px';cv.style.top=y+'px';W.appendChild(cv);}
for(let r=0;r<G.D-1;r++){for(const x of G.JCOL){
 const l=document.createElement('div');l.className='jlink';l.style.left=x+'px';
 l.style.top=(CY[r]+30)+'px';l.style.height=(CY[r+1]-CY[r]-60)+'px';W.appendChild(l);
 const d1=document.createElement('div');d1.className='jdot';d1.style.left=x+'px';d1.style.top=(CY[r]+30)+'px';W.appendChild(d1);
 const d2=document.createElement('div');d2.className='jdot';d2.style.left=x+'px';d2.style.top=(CY[r+1]-30)+'px';W.appendChild(d2);}}
const els={};
for(const i in KIND){const e=document.createElement('div');e.className='ion '+KIND[i];
 e.textContent=KIND[i]=='yb'?'':i;els[i]=e;W.appendChild(e);}
function jkey(s){return s.length==3?null:(Math.round(s[0]*10)/10)+'|'+s[1];}
function verify(fi){const f=FR[fi];const g={};
 for(const i in f.slots){const s=f.slots[i];if(s.length==2&&typeof s[0]=='number'){const k=jkey(s);(g[k]=g[k]||[]).push(i);}}
 for(const k in g){const x=parseFloat(k.split('|')[0]);
  if(CAP[x]!==undefined&&g[k].length>CAP[x])return 'capacity exceeded at x='+x;
  if(G.JCOL.some(j=>Math.abs(j-x)<1))return g[k]+' at rest on a junction column';}
 if(fi>0){const p=FR[fi-1];const share={};
  for(const fr of [p,f]){const h={};
   for(const i in fr.slots){const s=fr.slots[i];if(s.length==2&&typeof s[0]=='number'){const k=jkey(s);(h[k]=h[k]||[]).push(i);}}
   for(const k in h)for(const a of h[k])for(const b of h[k])if(a<b)share[a+'|'+b]=1;}
  for(let r=0;r<G.D;r++){
   const seq=fr=>Object.entries(fr.slots).filter(([i,s])=>s.length==2&&typeof s[0]=='number'&&s[1]==r)
     .sort((a,b)=>a[1][0]-b[1][0]).map(([i])=>i);
   const sa=seq(p),sb=seq(f),rk={};sb.forEach((i,k)=>rk[i]=k);
   const cm=sa.filter(i=>rk[i]!==undefined);
   for(let a=0;a<cm.length;a++)for(let b=a+1;b<cm.length;b++)
    if(rk[cm[a]]>rk[cm[b]]&&!share[[cm[a],cm[b]].sort().join('|')])
     return cm[a]+' passed '+cm[b]+' with no shared well (row '+r+')';}}
 return null;}
let cur=0,timer=null;
const cap=document.getElementById('cap'),ver=document.getElementById('verify'),
 sl=document.getElementById('sl'),ctr=document.getElementById('ctr');
sl.max=FR.length-1;
function show(k){cur=Math.max(0,Math.min(FR.length-1,k));const f=FR[cur];
 for(const i in f.pos){els[i].style.left=f.pos[i][0]+'px';els[i].style.top=f.pos[i][1]+'px';
  els[i].classList.toggle('hi',f.hi.includes(i));}
 cap.innerHTML=f.cap+(f.badge?'<span id="badge">'+f.badge+'</span>':'');
 const v=verify(cur);
 ver.innerHTML=v?'<span class="bad">&#9888; '+v+'</span>':'<span class="ok">&#10003; frame verified: capacities ok, junction columns clear, no free passing</span>';
 sl.value=cur;ctr.textContent=(cur+1)+' / '+FR.length;}
function play(){if(timer){clearInterval(timer);timer=null;document.getElementById('pl').textContent='Play';return;}
 document.getElementById('pl').textContent='Pause';
 timer=setInterval(()=>{if(cur>=FR.length-1){play();}else show(cur+1);},+document.getElementById('sp').value);}
document.getElementById('b0').onclick=()=>show(0);
document.getElementById('be').onclick=()=>show(FR.length-1);
document.getElementById('bp').onclick=()=>show(cur-1);
document.getElementById('bn').onclick=()=>show(cur+1);
document.getElementById('pl').onclick=play;
sl.oninput=e=>show(+e.target.value);
document.addEventListener('keydown',e=>{if(e.key=='ArrowRight')show(cur+1);if(e.key=='ArrowLeft')show(cur-1);if(e.key==' '){e.preventDefault();play();}});
show(0);
</script></body></html>"""
html = html.replace("__GEOM__", json.dumps(G)).replace("__FRAMES__", json.dumps(FR)).replace("__IONS__", json.dumps(IONS)).replace("__CAP__", json.dumps({str(k):v for k,v in CAPMAP.items()}).replace('"',''))
title = {"round":"Two local syndrome-extraction rounds, d = 3, replayed on the chapter-4 cell geometry",
         "merge":"One two-round remote-surgery merge, d = 3, replayed on the chapter-4 cell geometry"}[MODE]
html = html.replace("__TITLE__", title)
open(os.path.join(HERE, f"qec_chip_sim_d3_{MODE}.html"),"w").write(html)
print("written", len(html), "bytes")
