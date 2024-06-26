#   Copyright 2024 Muhammad Luckman Qasim,  Laleh Alisaraie
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

from jinja2 import Environment, FileSystemLoader
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.DSSP import DSSP
import math
import pdfkit
import base64
import cairosvg
import PIL
import imgkit
import pdf2image
import requests

PIL.Image.MAX_IMAGE_PIXELS = None


class VisualMap:
    '''
    Computes a Visualization mapping, residues with their secondary structure types.

    '''
    COLORS = {
        'H_COLOR': '#0000ff',
        'E_COLOR': '#ff0000',
        'E_A_COLOR': '#ff0000',
        'B_COLOR': '#228B22',
        'B_A_COLOR': '#228B22',
        'I_COLOR': '#ff00ff',
        'T_COLOR': '#ffff00',
        'S_COLOR': '#ffa500',
        'G_COLOR': '#ff0000',
        '-_COLOR': '#b7b7b7',
        'P_COLOR': '#b7b7b7'
    }
    ICONS = {
        'H': '<?xml version="1.0" encoding="UTF-8" standalone="no"?> <svg class="icon" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:cc="http://creativecommons.org/ns#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:svg="http://www.w3.org/2000/svg" xmlns="http://www.w3.org/2000/svg" version="1.1" id="svg4485" width="45" height="70" viewBox="0 0 45 70"> <defs id="defs4489" /> <path style="fill:{color};fill-opacity:1;stroke-width:0.64925498" d="M 37.848252,50.651288 C 35.158012,49.282652 31.538905,42.925169 29.054794,32.294924 27.52851,25.763522 25.383545,19.931385 24.507687,19.931385 c -1.039297,0 -5.93362,9.357192 -9.044078,17.290834 -1.57809,4.025146 -3.060396,6.39719 -4.368307,8.515175 C 7.676557,55.467756 3.289281e-6,55.539 3.289281e-6,55.539 L 0,41.045416 c -1.7361604e-6,-7.650057 0.00185329,-12.336357 0.00185329,-12.336357 0,0 1.65548821,3.498058 1.82624851,4.412888 0.2747429,1.471912 0.4181418,1.538248 1.2459953,0.576362 1.283722,-1.491553 2.4111662,-5.511281 2.4111662,-8.596657 0,-2.997608 2.3200041,-8.499697 4.8332857,-11.462494 1.354381,-1.596627 2.527089,-2.136718 5.440326,-2.505542 2.042576,-0.258602 5.513427,-0.805607 7.712994,-1.2155771 5.689088,-1.060357 8.31784,-0.8915145 10.698676,0.6871611 2.231106,1.479413 6.103048,8.000863 8.356539,14.074809 0.733694,1.977751 2.474771,4.029071 2.474771,4.029071 0,0 -0.02302,4.772834 -0.02302,11.158552 0,6.693024 0.02117,15.671387 0.02117,15.671387 0,0 -4.821633,-3.250321 -7.151553,-4.887731 z M 27.853519,19.084357 c 0.229909,-1.658974 -1.060782,-3.61579 -2.764683,-4.191528 -1.311304,-0.443084 -1.466021,-0.319774 -1.152446,0.918543 0.524251,2.070273 2.506403,4.963475 3.183999,4.64746 0.319072,-0.148808 0.648988,-0.767321 0.73313,-1.374475 z" id="path4510-8" /> <g id="g4530-3" transform="matrix(0.54870255,0,0,0.76848153,0.00126643,2.1515045)" style="opacity:0.456"> <path style="fill:#b9b9b9;fill-opacity:1" d="M 67.193339,62.127255 C 63.784542,60.483775 56.741729,50.710524 52.613775,39.456894 50.587618,33.933185 48.022545,27.921838 46.913613,26.098345 44.954554,22.876931 44.841237,22.833703 42.918065,24.57415 40.340859,26.90649 31.226244,40.660552 26.561435,49.256498 23.759798,54.419137 21.971078,56.484269 18.978906,58.010761 15.413576,59.829657 4.2955425,62.949206 1.25,62.985219 0.30915652,62.996344 0.00378724,59.599278 0.01531277,49.25 0.02373477,41.6875 0.00106954,34.558481 0.00106954,34.558481 c 0,0 1.58499126,2.750871 2.39664126,3.85113 0.3254234,0.441138 0.8108459,1.420605 1.0592729,2.035237 0.2334284,0.577525 0.3461697,1.449022 0.9954709,1.375557 0.9701687,-0.10977 2.484466,-2.860058 3.6175842,-4.714339 1.0263826,-1.679614 1.4004592,-7.38002 2.3217892,-9.831136 0.830557,-2.209622 1.635078,-4.289431 2.443584,-5.440184 4.296929,-6.115856 8.114404,-9.207842 11.389868,-9.696629 11.386772,-1.699212 27.933744,-3.0140136 30.950443,-3.0129532 3.710766,0.0013 9.906849,3.0114642 12.136544,5.7433582 4.210524,3.39774 11.637159,17.390392 13.821986,20.352449 C 82.052312,35.653348 82,39.666667 82,49 L 81.72727,62.774334 75.25,63.93232 C 71.5375,63.8951 68.4361,62.726427 67.193339,62.127255 Z" id="path4536-8" /> <path style="fill:#8c8c8c;fill-opacity:1" d="M 69,62.974804 C 63.368195,61.654591 57.985121,53.601464 52.974157,39 49.931967,30.135353 46.232751,23 44.679198,23 42.805903,23 33.865474,35.160993 28.5054,45 25.958564,49.675 22.8066,54.466281 21.501037,55.64729 15.0252,63.683907 13.708063,67.920507 -0.00230205,69.47141 L 0,50.475 C 0.001263,40.052778 0.24682699,38 1.5,38 1.6227654,38 2.8791855,39.231902 2.9320418,39.530078 3.446868,42.434338 8.3492855,51 10.196468,51 c 2.027976,0 2.388241,-2.015291 0.667818,-3.735714 -0.847619,-0.847619 -0.694054,-1.577375 0.605376,-2.876805 1.778307,-1.778307 1.77441,-2.597755 -0.0758,-15.941264 -0.50396,-3.634501 -0.140495,-4.587152 3.514122,-9.210589 5.083443,-6.431038 6.021303,-6.923385 14.616588,-7.673243 3.836487,-0.334698 11.306316,-1.113708 16.59962,-1.7311348 8.586429,-1.001546 10.171991,-0.9256065 14.704508,0.7042638 4.252785,1.52928 5.902182,2.92863 10.125811,8.590747 C 73.729529,22.846401 76,26.187106 76,26.550051 c 0,0.362945 1.35,2.66022 3,5.105056 2.961073,4.387478 3,4.618166 3,17.778394 0,12.981733 -1.012224,16.848237 -3.202904,17.435318 C 75.927588,67.637817 73.423257,64.011708 69,62.974804 Z M 54.75,19.814213 C 56.5375,17.04208 58,14.410582 58,13.966441 c 0,-1.344866 -5.778372,-2.176102 -10.271726,-1.477618 -8.231227,1.27953 -8.150147,1.151002 -3.987583,6.321097 5.496503,6.826911 5.6564,6.980604 6.754045,6.492009 C 51.047631,25.055818 52.9625,22.586346 54.75,19.814213 Z" id="path4534-5" /> <path style="fill:#606060;fill-opacity:1" d="M 66.5,61.792183 C 65.4,61.161334 63.955514,59.937518 63.290031,59.072592 61.476307,56.715302 55.99412,46.346346 56.040038,45.360023 56.062058,44.88701 56.722484,45.682254 57.50765,47.127232 59.398009,50.606147 62.86283,51.154873 72,49.422398 82.698986,47.393791 81.994297,46.968977 82,55.5 l 0.0093,13.971434 c 0,0 -9.080445,-4.406368 -15.50934,-7.679251 z M 0,54.5 c 0,-6.483035 0.23729179,-7.497564 1.75,-7.482033 0.9836017,0.0101 3.1673731,2.185748 4.9864409,4.967895 2.3373646,3.574848 4.1430511,5.18514 6.5000001,5.796633 3.18861,0.82726 3.310922,5.543595 -3.240599,7.190289 C -0.21120104,67.538276 0,63.125012 0,54.5 Z M 31.350286,37.088327 c 0.725447,-3.198242 -0.57457,-5.480463 -6.985964,-12.26409 -5.813986,-6.151533 -6.33122,-7.639232 -3.535581,-10.169249 2.145144,-1.941329 6.793284,-2.121829 12.555239,-0.487552 C 37.631545,15.372182 44,19.446514 44,20.959223 c 0,0.471592 -1.978272,3.285977 -4.39616,6.254188 -2.417888,2.968211 -5.387156,6.946962 -6.598374,8.841668 -1.784393,2.791326 -2.098426,2.987361 -1.65518,1.033248 z" id="path4532-0" /> </g> </svg>',
        'B': '<?xml version="1.0" encoding="UTF-8" standalone="no"?> <svg class="icon" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:cc="http://creativecommons.org/ns#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:svg="http://www.w3.org/2000/svg" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" version="1.1" id="svg4495" width="45" height="70" viewBox="0 0 45 70" sodipodi:docname="BetaStrand7.svg" inkscape:version="0.92.1 r15371"> <sodipodi:namedview pagecolor="#ffffff" bordercolor="#666666" borderopacity="1" objecttolerance="10" gridtolerance="10" guidetolerance="10" inkscape:pageopacity="0" inkscape:pageshadow="2" inkscape:window-width="1366" inkscape:window-height="705" id="namedview4578" showgrid="false" inkscape:zoom="3.55" inkscape:cx="-5.7746479" inkscape:cy="12.888617" inkscape:window-x="-8" inkscape:window-y="-8" inkscape:window-maximized="1" inkscape:current-layer="svg4495" /> <defs id="defs4499"> <linearGradient inkscape:collect="always" id="linearGradient5127"> <stop style="stop-color:{color};stop-opacity:1;" offset="0" id="stop5123" /> <stop style="stop-color:{color};stop-opacity:0;" offset="1" id="stop5125" /> </linearGradient> <linearGradient inkscape:collect="always" xlink:href="#linearGradient5127" id="linearGradient5129" x1="33.239433" y1="47.042252" x2="32.676056" y2="28.450703" gradientUnits="userSpaceOnUse" gradientTransform="matrix(0.56226913,0,0,1.427842,9.2818426e-6,-14.904152)" /> </defs> <path style="fill:url(#linearGradient5129);fill-opacity:1;stroke-width:0.11613007" d="m 44.981539,55.539131 c 0,0 -29.982112,0.458457 -44.97182625,3.73e-4 L 0.00185154,28.709003 44.983394,28.709377 Z" id="path5048" inkscape:connector-curvature="0" /> <g id="g5102" transform="matrix(0.11590555,0,0,0.11586807,-8.4137034,16.569976)" style="opacity:0.22300002"> <path style="fill:#aaaaaa;fill-opacity:1" d="m 460.67892,336.32373 -388.087704,0.003 0.01578,-231.56081 388.087924,0.003 z" id="path5108" inkscape:connector-curvature="0" /> <path style="fill:#000000;fill-opacity:1" d="m 460.67892,336.32373 -388.004089,0.003 0.07584,-53.35919 388.087709,-0.003 z" id="path5106" inkscape:connector-curvature="0" /> </g> </svg>',
        'B_A': '<?xml version="1.0" encoding="UTF-8" standalone="no"?> <svg class="icon" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:cc="http://creativecommons.org/ns#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:svg="http://www.w3.org/2000/svg" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" id="svg4495" width="45" height="70" viewBox="0 0 45 70"> <defs id="defs4499"> <linearGradient id="linearGradient5160"> <stop id="stop5156" offset="0" style="stop-color:{color};stop-opacity:1;" /> <stop id="stop5158" offset="1" style="stop-color:{color};stop-opacity:0;" /> </linearGradient> <linearGradient gradientTransform="matrix(0.56238502,0,0,1,0.01172915,9.9905989)" gradientUnits="userSpaceOnUse" y2="12" x2="32" y1="43" x1="32.75" id="linearGradient5162" xlink:href="#linearGradient5160" /> </defs> <path style="fill:url(#linearGradient5162);fill-opacity:1;stroke-width:0.08428686" d="M 24.21077,56.377974 C 16.137025,56.607263 8.0823627,55.977806 0.01171583,55.529617 L 0.0135707,28.699678 24.099421,29.360034 c 0.05979,-3.189974 0.07959,-10.561962 0.211293,-13.619267 7.305529,7.156614 14.229387,16.708881 20.620917,23.25133 0.06979,0.122159 -0.06113,0.662598 -0.774968,3.199786 -0.472402,1.679013 -0.890084,3.126005 -0.92816,3.215515 -6.718947,7.92766 -11.931461,16.633418 -18.843835,24.583201 -0.162223,-3.012579 -0.129919,-9.69331 -0.173898,-13.612625 z" id="path5048" /> <g id="g5102" transform="matrix(0.06240497,0,0,0.11336505,-4.4760318,17.493762)" style="opacity:0.21300001"> <path style="fill:#aaaaaa;fill-opacity:1" d="m 459.68774,343 c -129.07829,2.4063 -257.95058,-3.3447 -386.975186,-7.45438 l 0.108524,-42.63695 -1.065796,-192.5606 C 200.35448,107.45595 330.0348,106 458.7027,106 c 0.44776,-28.412551 2.26428,-120.629893 2.26428,-120.629893 0,0 234.84679,154.837803 314.53302,195.699413 9.59744,4.91981 17.13759,9.39809 17.31599,10.28438 0.17379,0.86344 -3.48404,8.89953 -8.1285,17.85799 -4.64447,8.95846 -11.03444,21.37586 -14.19994,27.59423 l -5.75545,11.30611 -17.07323,11.19389 C 641.9609,328.06111 559.66546,384.1368 462.50916,447 c -2.64853,-34.05546 -1.54843,-72.71659 -2.82142,-104 z" id="path5108" /> <path style="fill:#000000;fill-opacity:1" d="M 460.67903,343 283.58951,342.9847 C 186.19028,342.9767 71.725559,337.01659 71.725559,337.01659 L 72.021787,282.479 c 114.649933,-1.01054 385.286843,1.98867 385.286843,1.98867 z" id="path5106" /> <path style="fill:#000000;fill-opacity:1" d="m 463.00406,437.90499 c 0.002,-4.07275 0.29488,-13.10725 0.65035,-20.07668 l 0.6463,-12.67169 13.59964,-8.93971 C 588.70202,325.60067 677.89292,262.67778 791,193.2842 c 0.97687,-0.30795 1.11106,-0.1784 0.38478,0.37147 -0.61337,0.46438 -1.46152,1.96933 -1.88478,3.34433 -0.42326,1.375 -1.04995,2.725 -1.39266,3 -0.54593,0.43809 -18.93998,36.09649 -23.06707,44.71743 -0.84715,1.76959 -1.80983,3.23209 -2.13929,3.25 -103.01856,69.64434 -299.85716,215.19329 -299.85716,215.19329 0,0 -0.0437,-18.14513 -0.0398,-25.25573 z" id="path5104" /> </g> </svg>',
        'E': '<?xml version="1.0" encoding="UTF-8" standalone="no"?> <svg class="icon" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:cc="http://creativecommons.org/ns#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:svg="http://www.w3.org/2000/svg" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" version="1.1" id="svg4495" width="45" height="70" viewBox="0 0 45 70" sodipodi:docname="BetaStrand7.svg" inkscape:version="0.92.1 r15371"> <sodipodi:namedview pagecolor="#ffffff" bordercolor="#666666" borderopacity="1" objecttolerance="10" gridtolerance="10" guidetolerance="10" inkscape:pageopacity="0" inkscape:pageshadow="2" inkscape:window-width="1366" inkscape:window-height="705" id="namedview4578" showgrid="false" inkscape:zoom="3.55" inkscape:cx="-5.7746479" inkscape:cy="12.888617" inkscape:window-x="-8" inkscape:window-y="-8" inkscape:window-maximized="1" inkscape:current-layer="svg4495" /> <defs id="defs4499"> <linearGradient inkscape:collect="always" id="linearGradient5127"> <stop style="stop-color:{color};stop-opacity:1;" offset="0" id="stop5123" /> <stop style="stop-color:{color};stop-opacity:0;" offset="1" id="stop5125" /> </linearGradient> <linearGradient inkscape:collect="always" xlink:href="#linearGradient5127" id="linearGradient5129" x1="33.239433" y1="47.042252" x2="32.676056" y2="28.450703" gradientUnits="userSpaceOnUse" gradientTransform="matrix(0.56226913,0,0,1.427842,9.2818426e-6,-14.904152)" /> </defs> <path style="fill:url(#linearGradient5129);fill-opacity:1;stroke-width:0.11613007" d="m 44.981539,55.539131 c 0,0 -29.982112,0.458457 -44.97182625,3.73e-4 L 0.00185154,28.709003 44.983394,28.709377 Z" id="path5048" inkscape:connector-curvature="0" /> <g id="g5102" transform="matrix(0.11590555,0,0,0.11586807,-8.4137034,16.569976)" style="opacity:0.22300002"> <path style="fill:#aaaaaa;fill-opacity:1" d="m 460.67892,336.32373 -388.087704,0.003 0.01578,-231.56081 388.087924,0.003 z" id="path5108" inkscape:connector-curvature="0" /> <path style="fill:#000000;fill-opacity:1" d="m 460.67892,336.32373 -388.004089,0.003 0.07584,-53.35919 388.087709,-0.003 z" id="path5106" inkscape:connector-curvature="0" /> </g> </svg>',
        'E_A': '<?xml version="1.0" encoding="UTF-8" standalone="no"?> <svg class="icon" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:cc="http://creativecommons.org/ns#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:svg="http://www.w3.org/2000/svg" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" id="svg4495" width="45" height="70" viewBox="0 0 45 70"> <defs id="defs4499"> <linearGradient id="linearGradient5160"> <stop id="stop5156" offset="0" style="stop-color:{color};stop-opacity:1;" /> <stop id="stop5158" offset="1" style="stop-color:{color};stop-opacity:0;" /> </linearGradient> <linearGradient gradientTransform="matrix(0.56238502,0,0,1,0.01172915,9.9905989)" gradientUnits="userSpaceOnUse" y2="12" x2="32" y1="43" x1="32.75" id="linearGradient5162" xlink:href="#linearGradient5160" /> </defs> <path style="fill:url(#linearGradient5162);fill-opacity:1;stroke-width:0.08428686" d="M 24.21077,56.377974 C 16.137025,56.607263 8.0823627,55.977806 0.01171583,55.529617 L 0.0135707,28.699678 24.099421,29.360034 c 0.05979,-3.189974 0.07959,-10.561962 0.211293,-13.619267 7.305529,7.156614 14.229387,16.708881 20.620917,23.25133 0.06979,0.122159 -0.06113,0.662598 -0.774968,3.199786 -0.472402,1.679013 -0.890084,3.126005 -0.92816,3.215515 -6.718947,7.92766 -11.931461,16.633418 -18.843835,24.583201 -0.162223,-3.012579 -0.129919,-9.69331 -0.173898,-13.612625 z" id="path5048" /> <g id="g5102" transform="matrix(0.06240497,0,0,0.11336505,-4.4760318,17.493762)" style="opacity:0.21300001"> <path style="fill:#aaaaaa;fill-opacity:1" d="m 459.68774,343 c -129.07829,2.4063 -257.95058,-3.3447 -386.975186,-7.45438 l 0.108524,-42.63695 -1.065796,-192.5606 C 200.35448,107.45595 330.0348,106 458.7027,106 c 0.44776,-28.412551 2.26428,-120.629893 2.26428,-120.629893 0,0 234.84679,154.837803 314.53302,195.699413 9.59744,4.91981 17.13759,9.39809 17.31599,10.28438 0.17379,0.86344 -3.48404,8.89953 -8.1285,17.85799 -4.64447,8.95846 -11.03444,21.37586 -14.19994,27.59423 l -5.75545,11.30611 -17.07323,11.19389 C 641.9609,328.06111 559.66546,384.1368 462.50916,447 c -2.64853,-34.05546 -1.54843,-72.71659 -2.82142,-104 z" id="path5108" /> <path style="fill:#000000;fill-opacity:1" d="M 460.67903,343 283.58951,342.9847 C 186.19028,342.9767 71.725559,337.01659 71.725559,337.01659 L 72.021787,282.479 c 114.649933,-1.01054 385.286843,1.98867 385.286843,1.98867 z" id="path5106" /> <path style="fill:#000000;fill-opacity:1" d="m 463.00406,437.90499 c 0.002,-4.07275 0.29488,-13.10725 0.65035,-20.07668 l 0.6463,-12.67169 13.59964,-8.93971 C 588.70202,325.60067 677.89292,262.67778 791,193.2842 c 0.97687,-0.30795 1.11106,-0.1784 0.38478,0.37147 -0.61337,0.46438 -1.46152,1.96933 -1.88478,3.34433 -0.42326,1.375 -1.04995,2.725 -1.39266,3 -0.54593,0.43809 -18.93998,36.09649 -23.06707,44.71743 -0.84715,1.76959 -1.80983,3.23209 -2.13929,3.25 -103.01856,69.64434 -299.85716,215.19329 -299.85716,215.19329 0,0 -0.0437,-18.14513 -0.0398,-25.25573 z" id="path5104" /> </g> </svg>',
        'I': '<?xml version="1.0" encoding="UTF-8" standalone="no"?> <svg class="icon" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:cc="http://creativecommons.org/ns#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:svg="http://www.w3.org/2000/svg" xmlns="http://www.w3.org/2000/svg" version="1.1" id="svg4485" width="45" height="70" viewBox="0 0 45 70"> <defs id="defs4489" /> <path style="fill:{color};fill-opacity:1;stroke-width:0.64925498" d="M 37.848252,50.651288 C 35.158012,49.282652 31.538905,42.925169 29.054794,32.294924 27.52851,25.763522 25.383545,19.931385 24.507687,19.931385 c -1.039297,0 -5.93362,9.357192 -9.044078,17.290834 -1.57809,4.025146 -3.060396,6.39719 -4.368307,8.515175 C 7.676557,55.467756 3.289281e-6,55.539 3.289281e-6,55.539 L 0,41.045416 c -1.7361604e-6,-7.650057 0.00185329,-12.336357 0.00185329,-12.336357 0,0 1.65548821,3.498058 1.82624851,4.412888 0.2747429,1.471912 0.4181418,1.538248 1.2459953,0.576362 1.283722,-1.491553 2.4111662,-5.511281 2.4111662,-8.596657 0,-2.997608 2.3200041,-8.499697 4.8332857,-11.462494 1.354381,-1.596627 2.527089,-2.136718 5.440326,-2.505542 2.042576,-0.258602 5.513427,-0.805607 7.712994,-1.2155771 5.689088,-1.060357 8.31784,-0.8915145 10.698676,0.6871611 2.231106,1.479413 6.103048,8.000863 8.356539,14.074809 0.733694,1.977751 2.474771,4.029071 2.474771,4.029071 0,0 -0.02302,4.772834 -0.02302,11.158552 0,6.693024 0.02117,15.671387 0.02117,15.671387 0,0 -4.821633,-3.250321 -7.151553,-4.887731 z M 27.853519,19.084357 c 0.229909,-1.658974 -1.060782,-3.61579 -2.764683,-4.191528 -1.311304,-0.443084 -1.466021,-0.319774 -1.152446,0.918543 0.524251,2.070273 2.506403,4.963475 3.183999,4.64746 0.319072,-0.148808 0.648988,-0.767321 0.73313,-1.374475 z" id="path4510-8" /> <g id="g4530-3" transform="matrix(0.54870255,0,0,0.76848153,0.00126643,2.1515045)" style="opacity:0.456"> <path style="fill:#b9b9b9;fill-opacity:1" d="M 67.193339,62.127255 C 63.784542,60.483775 56.741729,50.710524 52.613775,39.456894 50.587618,33.933185 48.022545,27.921838 46.913613,26.098345 44.954554,22.876931 44.841237,22.833703 42.918065,24.57415 40.340859,26.90649 31.226244,40.660552 26.561435,49.256498 23.759798,54.419137 21.971078,56.484269 18.978906,58.010761 15.413576,59.829657 4.2955425,62.949206 1.25,62.985219 0.30915652,62.996344 0.00378724,59.599278 0.01531277,49.25 0.02373477,41.6875 0.00106954,34.558481 0.00106954,34.558481 c 0,0 1.58499126,2.750871 2.39664126,3.85113 0.3254234,0.441138 0.8108459,1.420605 1.0592729,2.035237 0.2334284,0.577525 0.3461697,1.449022 0.9954709,1.375557 0.9701687,-0.10977 2.484466,-2.860058 3.6175842,-4.714339 1.0263826,-1.679614 1.4004592,-7.38002 2.3217892,-9.831136 0.830557,-2.209622 1.635078,-4.289431 2.443584,-5.440184 4.296929,-6.115856 8.114404,-9.207842 11.389868,-9.696629 11.386772,-1.699212 27.933744,-3.0140136 30.950443,-3.0129532 3.710766,0.0013 9.906849,3.0114642 12.136544,5.7433582 4.210524,3.39774 11.637159,17.390392 13.821986,20.352449 C 82.052312,35.653348 82,39.666667 82,49 L 81.72727,62.774334 75.25,63.93232 C 71.5375,63.8951 68.4361,62.726427 67.193339,62.127255 Z" id="path4536-8" /> <path style="fill:#8c8c8c;fill-opacity:1" d="M 69,62.974804 C 63.368195,61.654591 57.985121,53.601464 52.974157,39 49.931967,30.135353 46.232751,23 44.679198,23 42.805903,23 33.865474,35.160993 28.5054,45 25.958564,49.675 22.8066,54.466281 21.501037,55.64729 15.0252,63.683907 13.708063,67.920507 -0.00230205,69.47141 L 0,50.475 C 0.001263,40.052778 0.24682699,38 1.5,38 1.6227654,38 2.8791855,39.231902 2.9320418,39.530078 3.446868,42.434338 8.3492855,51 10.196468,51 c 2.027976,0 2.388241,-2.015291 0.667818,-3.735714 -0.847619,-0.847619 -0.694054,-1.577375 0.605376,-2.876805 1.778307,-1.778307 1.77441,-2.597755 -0.0758,-15.941264 -0.50396,-3.634501 -0.140495,-4.587152 3.514122,-9.210589 5.083443,-6.431038 6.021303,-6.923385 14.616588,-7.673243 3.836487,-0.334698 11.306316,-1.113708 16.59962,-1.7311348 8.586429,-1.001546 10.171991,-0.9256065 14.704508,0.7042638 4.252785,1.52928 5.902182,2.92863 10.125811,8.590747 C 73.729529,22.846401 76,26.187106 76,26.550051 c 0,0.362945 1.35,2.66022 3,5.105056 2.961073,4.387478 3,4.618166 3,17.778394 0,12.981733 -1.012224,16.848237 -3.202904,17.435318 C 75.927588,67.637817 73.423257,64.011708 69,62.974804 Z M 54.75,19.814213 C 56.5375,17.04208 58,14.410582 58,13.966441 c 0,-1.344866 -5.778372,-2.176102 -10.271726,-1.477618 -8.231227,1.27953 -8.150147,1.151002 -3.987583,6.321097 5.496503,6.826911 5.6564,6.980604 6.754045,6.492009 C 51.047631,25.055818 52.9625,22.586346 54.75,19.814213 Z" id="path4534-5" /> <path style="fill:#606060;fill-opacity:1" d="M 66.5,61.792183 C 65.4,61.161334 63.955514,59.937518 63.290031,59.072592 61.476307,56.715302 55.99412,46.346346 56.040038,45.360023 56.062058,44.88701 56.722484,45.682254 57.50765,47.127232 59.398009,50.606147 62.86283,51.154873 72,49.422398 82.698986,47.393791 81.994297,46.968977 82,55.5 l 0.0093,13.971434 c 0,0 -9.080445,-4.406368 -15.50934,-7.679251 z M 0,54.5 c 0,-6.483035 0.23729179,-7.497564 1.75,-7.482033 0.9836017,0.0101 3.1673731,2.185748 4.9864409,4.967895 2.3373646,3.574848 4.1430511,5.18514 6.5000001,5.796633 3.18861,0.82726 3.310922,5.543595 -3.240599,7.190289 C -0.21120104,67.538276 0,63.125012 0,54.5 Z M 31.350286,37.088327 c 0.725447,-3.198242 -0.57457,-5.480463 -6.985964,-12.26409 -5.813986,-6.151533 -6.33122,-7.639232 -3.535581,-10.169249 2.145144,-1.941329 6.793284,-2.121829 12.555239,-0.487552 C 37.631545,15.372182 44,19.446514 44,20.959223 c 0,0.471592 -1.978272,3.285977 -4.39616,6.254188 -2.417888,2.968211 -5.387156,6.946962 -6.598374,8.841668 -1.784393,2.791326 -2.098426,2.987361 -1.65518,1.033248 z" id="path4532-0" /> </g> </svg>',
        'T': '<?xml version="1.0" encoding="UTF-8" standalone="no"?> <svg class="icon" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:cc="http://creativecommons.org/ns#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:svg="http://www.w3.org/2000/svg" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 45 70" height="70" width="45" id="svg4495" version="1.1"> <defs id="defs4499"> <linearGradient id="linearGradient5127"> <stop id="stop5123" offset="0" style="stop-color:{color};stop-opacity:1;" /> <stop id="stop5125" offset="1" style="stop-color:{color};stop-opacity:0;" /> </linearGradient> <linearGradient gradientTransform="matrix(0.56226913,0,0,1.0563306,9.2818426e-6,-1.2654213)" gradientUnits="userSpaceOnUse" y2="28.450703" x2="32.676056" y1="47.042252" x1="33.239433" id="linearGradient5129" xlink:href="#linearGradient5127" /> <linearGradient xlink:href="#linearGradient5127" id="linearGradient5162" x1="32.75" y1="43" x2="32" y2="12" gradientUnits="userSpaceOnUse" gradientTransform="matrix(0.56238502,0,0,1,-44.988271,9.9905988)" /> </defs> <path id="path5048" d="m 44.981539,50.849176 c 0,0 -29.982112,0.339171 -44.97182625,2.76e-4 L 0.00185154,31.000002 44.983394,31.000278 Z" style="fill:url(#linearGradient5129);fill-opacity:1;stroke-width:0.09988598" /> <g style="opacity:0.22300002" transform="matrix(0.11590555,0,0,0.08572026,-8.4137034,22.019439)" id="g5102"> <path id="path5108" d="m 460.67892,336.32373 -388.087704,0.003 0.01578,-231.56081 388.087924,0.003 z" style="fill:#aaaaaa;fill-opacity:1" /> <path id="path5106" d="m 460.67892,336.32373 -388.004089,0.003 0.07584,-53.35919 388.087709,-0.003 z" style="fill:#000000;fill-opacity:1" /> </g> </svg>',
        'S': '<?xml version="1.0" encoding="UTF-8" standalone="no"?> <svg class="icon" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:cc="http://creativecommons.org/ns#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:svg="http://www.w3.org/2000/svg" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 45 70" height="70" width="45" id="svg4495" version="1.1"> <defs id="defs4499"> <linearGradient id="linearGradient5127"> <stop id="stop5123" offset="0" style="stop-color:{color};stop-opacity:1;" /> <stop id="stop5125" offset="1" style="stop-color:{color};stop-opacity:0;" /> </linearGradient> <linearGradient gradientTransform="matrix(0.56226913,0,0,1.0563306,9.2818426e-6,-1.2654213)" gradientUnits="userSpaceOnUse" y2="28.450703" x2="32.676056" y1="47.042252" x1="33.239433" id="linearGradient5129" xlink:href="#linearGradient5127" /> <linearGradient xlink:href="#linearGradient5127" id="linearGradient5162" x1="32.75" y1="43" x2="32" y2="12" gradientUnits="userSpaceOnUse" gradientTransform="matrix(0.56238502,0,0,1,-44.988271,9.9905988)" /> </defs> <path id="path5048" d="m 44.981539,50.849176 c 0,0 -29.982112,0.339171 -44.97182625,2.76e-4 L 0.00185154,31.000002 44.983394,31.000278 Z" style="fill:url(#linearGradient5129);fill-opacity:1;stroke-width:0.09988598" /> <g style="opacity:0.22300002" transform="matrix(0.11590555,0,0,0.08572026,-8.4137034,22.019439)" id="g5102"> <path id="path5108" d="m 460.67892,336.32373 -388.087704,0.003 0.01578,-231.56081 388.087924,0.003 z" style="fill:#aaaaaa;fill-opacity:1" /> <path id="path5106" d="m 460.67892,336.32373 -388.004089,0.003 0.07584,-53.35919 388.087709,-0.003 z" style="fill:#000000;fill-opacity:1" /> </g> </svg>',
        'G': '<?xml version="1.0" encoding="UTF-8" standalone="no"?> <svg class="icon" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:cc="http://creativecommons.org/ns#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:svg="http://www.w3.org/2000/svg" xmlns="http://www.w3.org/2000/svg" version="1.1" id="svg4485" width="45" height="70" viewBox="0 0 45 70"> <defs id="defs4489" /> <path style="fill:{color};fill-opacity:1;stroke-width:0.64925498" d="M 37.848252,50.651288 C 35.158012,49.282652 31.538905,42.925169 29.054794,32.294924 27.52851,25.763522 25.383545,19.931385 24.507687,19.931385 c -1.039297,0 -5.93362,9.357192 -9.044078,17.290834 -1.57809,4.025146 -3.060396,6.39719 -4.368307,8.515175 C 7.676557,55.467756 3.289281e-6,55.539 3.289281e-6,55.539 L 0,41.045416 c -1.7361604e-6,-7.650057 0.00185329,-12.336357 0.00185329,-12.336357 0,0 1.65548821,3.498058 1.82624851,4.412888 0.2747429,1.471912 0.4181418,1.538248 1.2459953,0.576362 1.283722,-1.491553 2.4111662,-5.511281 2.4111662,-8.596657 0,-2.997608 2.3200041,-8.499697 4.8332857,-11.462494 1.354381,-1.596627 2.527089,-2.136718 5.440326,-2.505542 2.042576,-0.258602 5.513427,-0.805607 7.712994,-1.2155771 5.689088,-1.060357 8.31784,-0.8915145 10.698676,0.6871611 2.231106,1.479413 6.103048,8.000863 8.356539,14.074809 0.733694,1.977751 2.474771,4.029071 2.474771,4.029071 0,0 -0.02302,4.772834 -0.02302,11.158552 0,6.693024 0.02117,15.671387 0.02117,15.671387 0,0 -4.821633,-3.250321 -7.151553,-4.887731 z M 27.853519,19.084357 c 0.229909,-1.658974 -1.060782,-3.61579 -2.764683,-4.191528 -1.311304,-0.443084 -1.466021,-0.319774 -1.152446,0.918543 0.524251,2.070273 2.506403,4.963475 3.183999,4.64746 0.319072,-0.148808 0.648988,-0.767321 0.73313,-1.374475 z" id="path4510-8" /> <g id="g4530-3" transform="matrix(0.54870255,0,0,0.76848153,0.00126643,2.1515045)" style="opacity:0.456"> <path style="fill:#b9b9b9;fill-opacity:1" d="M 67.193339,62.127255 C 63.784542,60.483775 56.741729,50.710524 52.613775,39.456894 50.587618,33.933185 48.022545,27.921838 46.913613,26.098345 44.954554,22.876931 44.841237,22.833703 42.918065,24.57415 40.340859,26.90649 31.226244,40.660552 26.561435,49.256498 23.759798,54.419137 21.971078,56.484269 18.978906,58.010761 15.413576,59.829657 4.2955425,62.949206 1.25,62.985219 0.30915652,62.996344 0.00378724,59.599278 0.01531277,49.25 0.02373477,41.6875 0.00106954,34.558481 0.00106954,34.558481 c 0,0 1.58499126,2.750871 2.39664126,3.85113 0.3254234,0.441138 0.8108459,1.420605 1.0592729,2.035237 0.2334284,0.577525 0.3461697,1.449022 0.9954709,1.375557 0.9701687,-0.10977 2.484466,-2.860058 3.6175842,-4.714339 1.0263826,-1.679614 1.4004592,-7.38002 2.3217892,-9.831136 0.830557,-2.209622 1.635078,-4.289431 2.443584,-5.440184 4.296929,-6.115856 8.114404,-9.207842 11.389868,-9.696629 11.386772,-1.699212 27.933744,-3.0140136 30.950443,-3.0129532 3.710766,0.0013 9.906849,3.0114642 12.136544,5.7433582 4.210524,3.39774 11.637159,17.390392 13.821986,20.352449 C 82.052312,35.653348 82,39.666667 82,49 L 81.72727,62.774334 75.25,63.93232 C 71.5375,63.8951 68.4361,62.726427 67.193339,62.127255 Z" id="path4536-8" /> <path style="fill:#8c8c8c;fill-opacity:1" d="M 69,62.974804 C 63.368195,61.654591 57.985121,53.601464 52.974157,39 49.931967,30.135353 46.232751,23 44.679198,23 42.805903,23 33.865474,35.160993 28.5054,45 25.958564,49.675 22.8066,54.466281 21.501037,55.64729 15.0252,63.683907 13.708063,67.920507 -0.00230205,69.47141 L 0,50.475 C 0.001263,40.052778 0.24682699,38 1.5,38 1.6227654,38 2.8791855,39.231902 2.9320418,39.530078 3.446868,42.434338 8.3492855,51 10.196468,51 c 2.027976,0 2.388241,-2.015291 0.667818,-3.735714 -0.847619,-0.847619 -0.694054,-1.577375 0.605376,-2.876805 1.778307,-1.778307 1.77441,-2.597755 -0.0758,-15.941264 -0.50396,-3.634501 -0.140495,-4.587152 3.514122,-9.210589 5.083443,-6.431038 6.021303,-6.923385 14.616588,-7.673243 3.836487,-0.334698 11.306316,-1.113708 16.59962,-1.7311348 8.586429,-1.001546 10.171991,-0.9256065 14.704508,0.7042638 4.252785,1.52928 5.902182,2.92863 10.125811,8.590747 C 73.729529,22.846401 76,26.187106 76,26.550051 c 0,0.362945 1.35,2.66022 3,5.105056 2.961073,4.387478 3,4.618166 3,17.778394 0,12.981733 -1.012224,16.848237 -3.202904,17.435318 C 75.927588,67.637817 73.423257,64.011708 69,62.974804 Z M 54.75,19.814213 C 56.5375,17.04208 58,14.410582 58,13.966441 c 0,-1.344866 -5.778372,-2.176102 -10.271726,-1.477618 -8.231227,1.27953 -8.150147,1.151002 -3.987583,6.321097 5.496503,6.826911 5.6564,6.980604 6.754045,6.492009 C 51.047631,25.055818 52.9625,22.586346 54.75,19.814213 Z" id="path4534-5" /> <path style="fill:#606060;fill-opacity:1" d="M 66.5,61.792183 C 65.4,61.161334 63.955514,59.937518 63.290031,59.072592 61.476307,56.715302 55.99412,46.346346 56.040038,45.360023 56.062058,44.88701 56.722484,45.682254 57.50765,47.127232 59.398009,50.606147 62.86283,51.154873 72,49.422398 82.698986,47.393791 81.994297,46.968977 82,55.5 l 0.0093,13.971434 c 0,0 -9.080445,-4.406368 -15.50934,-7.679251 z M 0,54.5 c 0,-6.483035 0.23729179,-7.497564 1.75,-7.482033 0.9836017,0.0101 3.1673731,2.185748 4.9864409,4.967895 2.3373646,3.574848 4.1430511,5.18514 6.5000001,5.796633 3.18861,0.82726 3.310922,5.543595 -3.240599,7.190289 C -0.21120104,67.538276 0,63.125012 0,54.5 Z M 31.350286,37.088327 c 0.725447,-3.198242 -0.57457,-5.480463 -6.985964,-12.26409 -5.813986,-6.151533 -6.33122,-7.639232 -3.535581,-10.169249 2.145144,-1.941329 6.793284,-2.121829 12.555239,-0.487552 C 37.631545,15.372182 44,19.446514 44,20.959223 c 0,0.471592 -1.978272,3.285977 -4.39616,6.254188 -2.417888,2.968211 -5.387156,6.946962 -6.598374,8.841668 -1.784393,2.791326 -2.098426,2.987361 -1.65518,1.033248 z" id="path4532-0" /> </g> </svg>',
        '-': '<?xml version="1.0" encoding="UTF-8" standalone="no"?> <svg class="icon" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:cc="http://creativecommons.org/ns#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:svg="http://www.w3.org/2000/svg" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 45 70" height="70" width="45" id="svg4495" version="1.1"> <defs id="defs4499"> <linearGradient id="linearGradient5127"> <stop id="stop5123" offset="0" style="stop-color:{color};stop-opacity:1;" /> <stop id="stop5125" offset="1" style="stop-color:{color};stop-opacity:0;" /> </linearGradient> <linearGradient gradientTransform="matrix(0.56226913,0,0,0.47534886,9.2818426e-6,23.623543)" gradientUnits="userSpaceOnUse" y2="28.450703" x2="32.676056" y1="47.042252" x1="33.239433" id="linearGradient5129" xlink:href="#linearGradient5127" /> <linearGradient xlink:href="#linearGradient5127" id="linearGradient5162" x1="32.75" y1="43" x2="32" y2="12" gradientUnits="userSpaceOnUse" gradientTransform="matrix(0.56238502,0,0,1,-44.988271,9.9905988)" /> <linearGradient xlink:href="#linearGradient5127" id="linearGradient5162-6" x1="32.75" y1="43" x2="32" y2="12" gradientUnits="userSpaceOnUse" gradientTransform="matrix(0.56238502,0,0,1,45.01173,9.9905986)" /> </defs> <path id="path5048" d="m 44.981539,47.075124 c 0,0 -29.982112,0.152632 -44.97182625,1.28e-4 L 0.00185154,38.142998 44.983394,38.143122 Z" style="fill:url(#linearGradient5129);fill-opacity:1;stroke-width:0.06700556" /> <g style="opacity:0.22300002" transform="matrix(0.11590555,0,0,0.03857415,-8.4137034,34.101743)" id="g5102"> <path id="path5108" d="m 460.67892,336.32373 -388.087704,0.003 0.01578,-231.56081 388.087924,0.003 z" style="fill:#aaaaaa;fill-opacity:1" /> <path id="path5106" d="m 460.67892,336.32373 -388.004089,0.003 0.07584,-53.35919 388.087709,-0.003 z" style="fill:#000000;fill-opacity:1" /> </g> </svg>',
        'P': '<?xml version="1.0" encoding="UTF-8" standalone="no"?> <svg class="icon" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:cc="http://creativecommons.org/ns#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:svg="http://www.w3.org/2000/svg" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 45 70" height="70" width="45" id="svg4495" version="1.1"> <defs id="defs4499"> <linearGradient id="linearGradient5127"> <stop id="stop5123" offset="0" style="stop-color:{color};stop-opacity:1;" /> <stop id="stop5125" offset="1" style="stop-color:{color};stop-opacity:0;" /> </linearGradient> <linearGradient gradientTransform="matrix(0.56226913,0,0,0.47534886,9.2818426e-6,23.623543)" gradientUnits="userSpaceOnUse" y2="28.450703" x2="32.676056" y1="47.042252" x1="33.239433" id="linearGradient5129" xlink:href="#linearGradient5127" /> <linearGradient xlink:href="#linearGradient5127" id="linearGradient5162" x1="32.75" y1="43" x2="32" y2="12" gradientUnits="userSpaceOnUse" gradientTransform="matrix(0.56238502,0,0,1,-44.988271,9.9905988)" /> <linearGradient xlink:href="#linearGradient5127" id="linearGradient5162-6" x1="32.75" y1="43" x2="32" y2="12" gradientUnits="userSpaceOnUse" gradientTransform="matrix(0.56238502,0,0,1,45.01173,9.9905986)" /> </defs> <path id="path5048" d="m 44.981539,47.075124 c 0,0 -29.982112,0.152632 -44.97182625,1.28e-4 L 0.00185154,38.142998 44.983394,38.143122 Z" style="fill:url(#linearGradient5129);fill-opacity:1;stroke-width:0.06700556" /> <g style="opacity:0.22300002" transform="matrix(0.11590555,0,0,0.03857415,-8.4137034,34.101743)" id="g5102"> <path id="path5108" d="m 460.67892,336.32373 -388.087704,0.003 0.01578,-231.56081 388.087924,0.003 z" style="fill:#aaaaaa;fill-opacity:1" /> <path id="path5106" d="m 460.67892,336.32373 -388.004089,0.003 0.07584,-53.35919 388.087709,-0.003 z" style="fill:#000000;fill-opacity:1" /> </g> </svg>'
    }
    
    def __init__(self, file_path: str, pdb_name: str = None, subtitle: str = None, scientific_name: str = None) -> None:
        '''
        
        Args:
            pdb_name (str): The PDB code of the file.
            file_path (str): The file path of the PDB file.

        '''
        self.structure_list = self._get_dssp_output(pdb_name, file_path)
        self.file_path = file_path
        self.pdb_name = pdb_name
        self.subtitle = subtitle
        self.scientific_name = scientific_name

    def _update_template(self, residues_per_line: int = 50) -> str:
        '''
        A private Method that updates the template and returns a string representing the updated template.

        '''
        environment = Environment(loader=FileSystemLoader("templates/"))
        template = environment.get_template("template.html.jinja")
        unitprot_data = self.get_uniprot_data(self.pdb_name)

        residue_tables = ''

        for chain_id, chain in self.structure_list.items():
            table = ''
            struc_annotation_row = ''
            structure_row = ''
            residue_row = ''
            last_count = 0
            curr_structure_start = 0
            index_count = {'H': 0, 'B': 0, 'E': 0, 'I': 0, 'G': 0, 'S': 0, 'T': 0}

            for i, res in enumerate(chain):
                # Get the structure icon
                if (i == (len(chain) - 1) or chain[i]['res_struc'] != chain[i+1]['res_struc']) and chain[i]['res_struc'] in ['B', 'E']:
                    bytestring = VisualMap.ICONS[f"{res['res_struc']}_A"].format(color=VisualMap.COLORS[f"{res['res_struc']}_A_COLOR"])
                    structure_icon = self._convert_to_png(bytestring)
                else:
                    bytestring = VisualMap.ICONS[res['res_struc']].format(color=VisualMap.COLORS[f"{res['res_struc']}_COLOR"])
                    structure_icon = self._convert_to_png(bytestring)

                # Create the table cells for the indexed secondary structure assignment
                if i == (len(chain) - 1) or chain[i]['res_struc'] != chain[i+1]['res_struc'] or (i + 1) % residues_per_line == 0:
                    col_span = i + 1 - curr_structure_start
                    curr_structure_start = i + 1
                    if chain[i]["res_struc"] in index_count.keys() and chain[i]['res_struc'] != chain[i+1]['res_struc']:
                        index_count[chain[i]['res_struc']] += 1
                        struc_annotation_row += f'<td id="annotation" colspan="{col_span}">{chain[i]["res_struc"]}{index_count[chain[i]["res_struc"]]}</td>' 
                    else:
                        struc_annotation_row += f'<td id="annotation" colspan="{col_span}"></td>'
                structure_row += f'<td id="icon_row">{structure_icon}</td>'
                residue_row += f'<td id="res_row">{res["res_name"]}</td>'

                if (i + 1) % residues_per_line == 0 or i == (len(chain) - 1):
                    count = i
                    while (count + 1) % residues_per_line != 0:
                        structure_row += f'<td></td>'
                        residue_row += f'<td></td>'
                        struc_annotation_row += f'<td id="annotation"></td>'
                        count += 1
                    table += f'<tr id="annotation_row"><td></td>{struc_annotation_row}<td></td></tr>'
                    table += f'<tr><td id="icon_row"></td>{structure_row}<td></td></tr>'
                    table += f'<tr><td id="count_element_start">{last_count + 1}</td>{residue_row}<td id="count_element_end">{i + 1}</td></tr>'
                    struc_annotation_row = ''
                    structure_row = ''
                    residue_row = ''
                    last_count = i + 1
            
            uniprot_id = self.get_uniprot_id_by_chain_id(unitprot_data, chain_id)
            residue_tables += f'<div id="chain">Chain {chain_id}:</div>'
            if uniprot_id:
                residue_tables += f'<div>Uniprot ID: {uniprot_id}</div>'
            residue_tables += f'<table id="res"><tbody id="res">{table}</tbody></table>'
            last_count = 0
            index_count = {'H': 0, 'B': 0}

        rcsb_data = self.get_rcsb_entry_data(self.pdb_name)
        if self.pdb_name:
            self.pdb_name = self.pdb_name.upper()        
        if self.subtitle == None and rcsb_data:
            self.subtitle = rcsb_data['struct']['title']
        if self.scientific_name == None:
            self.scientific_name = self.get_scientific_name(self.pdb_name, rcsb_data)

        content = template.render(
            pdb_name = self.pdb_name,
            pdb_title = self.subtitle,
            scientific_name = self.scientific_name,
            residue_tables = residue_tables,
            H = self._convert_to_png(VisualMap.ICONS['H'].format(color=VisualMap.COLORS['H_COLOR'])),
            B = self._convert_to_png(VisualMap.ICONS['B_A'].format(color=VisualMap.COLORS['B_A_COLOR'])),
            E = self._convert_to_png(VisualMap.ICONS['E_A'].format(color=VisualMap.COLORS['E_A_COLOR'])),
            G = self._convert_to_png(VisualMap.ICONS['G'].format(color=VisualMap.COLORS['G_COLOR'])),
            I = self._convert_to_png(VisualMap.ICONS['I'].format(color=VisualMap.COLORS['I_COLOR'])),
            T = self._convert_to_png(VisualMap.ICONS['T'].format(color=VisualMap.COLORS['T_COLOR'])),
            S = self._convert_to_png(VisualMap.ICONS['S'].format(color=VisualMap.COLORS['S_COLOR'])),
            P = self._convert_to_png(VisualMap.ICONS['P'].format(color=VisualMap.COLORS['P_COLOR']))
        )
            
        return content
    
    def _get_width_and_height(self, residues_per_line):
        '''
        Private method that returns the width and height of the resulting visualization in the form of a table. Eg. (output_width, output_height)

        '''
        output_width = (20 * residues_per_line + 48 + 100 + 50) / 2
        output_height = (len(self.structure_list) * 76 + sum(math.ceil(len(i) / residues_per_line) for i in self.structure_list.values()) * 96 + 100 + 77 + 79 + 400) / 2
        return (int(output_width), int(output_height))
    
    def _convert_to_png(self, bytestring):
        img = base64.b64encode(cairosvg.svg2png(bytestring=bytestring))
        structure_icon_png = f'<img class="icon" src="data:image/png;base64, {img.decode("utf-8")}" />'
        return structure_icon_png
    
    def _get_dssp_output(self, pdb_name, file_path):
        '''
        Private method that takes in the PDB name and path and returns the DSSP secondary structure assignments on the PDB file.

        Returns:
            A dictionary containing chain ID as the key, the values are a list of dictionaries. Each dictionary item in the list represents a residue, along with its secondary structure assignments.
            Example output: {"A": [{"chain": "A", "res_num": 1, "res_name": "A", "res_struc": "H"}, {}, {}]}

        '''
        if file_path.split('.')[-1].lower() == 'pdb':
            p = PDBParser(QUIET=True)
        elif file_path.split('.')[-1].lower() == 'cif':
            p = MMCIFParser(QUIET=True)
        else:
            return Exception('File type has to be either PDB or mmCIF')
        
        structure = p.get_structure(pdb_name, file_path)
        model = structure[0]
        dssp = DSSP(model, file_path)

        structure_out = {}
        for i in list(dssp.keys()):
            res_dict = {'chain': i[0][0], 'res_num': i[1][1], 'res_name': dssp[i][1], 'res_struc': dssp[i][2]}
            if i[0][0] not in structure_out:
                structure_out[i[0][0]] = []
            structure_out[i[0][0]].append(res_dict)

        return structure_out
    
    def get_uniprot_data(self, identifier):
        url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{identifier}"
        response = requests.get(url)
        
        if response.status_code == 200:
            json_data = response.json() 
            return json_data
        else:
            print(f"Error: Unable to fetch data (status code: {response.status_code})")
            return None
        
    def get_rcsb_entry_data(self, identifier):
        url = f"https://data.rcsb.org/rest/v1/core/entry/{identifier}"
        response = requests.get(url)
        
        if response.status_code == 200:
            json_data = response.json() 
            return json_data
        else:
            print(f"Error: Unable to fetch data (status code: {response.status_code})")
            return None
        
    def get_scientific_name(self, identifier, rcsb_data):
        if identifier and rcsb_data:
            entity_ids = rcsb_data['rcsb_entry_container_identifiers']['polymer_entity_ids']

            for entity_id in entity_ids:
                url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{identifier}/{entity_id}"
                response = requests.get(url)
                
                if response.status_code == 200:
                    json_data = response.json() 
                else:
                    print(f"Error: Unable to fetch data (status code: {response.status_code})")
                    continue
                
                if 'rcsb_entity_source_organism' in json_data:
                    return json_data['rcsb_entity_source_organism'][0].get('scientific_name', None)
                else:
                    continue
        return None
    
    def get_uniprot_id_by_chain_id(self, data, target_chain_id):
        if data and target_chain_id:
            for structure_id, structure_data in data.items():
                for uniprot_id, uniprot_data in structure_data['UniProt'].items():
                    for mapping in uniprot_data['mappings']:
                        if mapping['chain_id'] == target_chain_id:
                            return uniprot_id
        return None

    def generate_visual(self, residues_per_line: int = 50, output_image_name: str = '', dpi: int = 100, pdf: bool = False) -> None:
        '''
        Generates a visualization with the residues mapped to their secondar structure types.

        Args:
            residues_per_line (int): An integer respresnting the number of residues per line, in the visualization, the default values is 60.
            output_image_name (string): Representing the output image path, valid types: "jpg, png, and gif".
            pdf (bool): A boolean value indicating if a PDF should be generated too. The default value is False.
            color: TO BE ADDED

        Raises:
            TypeError: If the any of the parameters are not of the correct type
            ValueError: If any of the parameters are outside the range

        '''
        if not isinstance(residues_per_line, int):
            raise TypeError('The residues_per_line parameter must be an integer')
        if not isinstance(output_image_name, str):
            raise TypeError('The ouput_image parameter must be a string')
        if not isinstance(pdf, bool):
            raise TypeError('The pdf parameter must be a boolean value (True/False)')

        if output_image_name == '':
            output_image_name = f'{self.pdb_name}.png'
        if '.' not in output_image_name or output_image_name.split('.')[-1].lower() not in ['jpg', 'png']:
            raise ValueError('The output image name must end with one of the following extensions: "JPG", "PNG"')
        
        updated_template = self._update_template(residues_per_line=residues_per_line)
        output_width, output_height = self._get_width_and_height(residues_per_line=residues_per_line)
        output_width, output_height = self._get_width_and_height(residues_per_line=residues_per_line)

        options = {'page-height': f'{output_height}px',
                   'page-width': f'{output_width}px',
                   'margin-top': '0',
                   'margin-bottom': '0',
                   'margin-left': '0',
                   'margin-right': '0',
                   'enable-local-file-access': '',
                   'print-media-type': True,
                   'disable-smart-shrinking': True}
        if dpi == 100:
            imgkit.from_string(updated_template, f'{output_image_name}', css='templates/output_styles.css', options={'quality': 100})
        else:
            pdf_bytes = pdfkit.from_string(updated_template, css='templates/output_styles.css', options=options)
            pages = pdf2image.convert_from_bytes(pdf_bytes, dpi=dpi)
            if len(pages) == 1:
                pages[0].save(output_image_name)
            else:
                for count, page in enumerate(pages):
                    page.save(f'{output_image_name}')

        if pdf:
            pdfkit.from_string(updated_template, f'{self.pdb_name}.pdf', css='templates/output_styles.css', options=options)
