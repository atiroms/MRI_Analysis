sudo docker run --rm -it -v /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC1_T1/CSUB-{subject}C-01/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic_T1.py -s 00003 00004 -ss 001 -c dcm2niix -b --overwrite


sudo docker run --rm -it -v /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC1_*/CSUB-{subject}C-01/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00003 00004 -ss 001 -c dcm2niix -b --overwrite


sudo docker run --rm -it -v /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC2_*/CSUB-{subject}C-02/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00003 00004 -ss 002 -c dcm2niix -b --overwrite


##
sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS/test2:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC1_T1/CSUB-{subject}C-01/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00014 00019 -ss 01 -c dcm2niix -b --overwrite

sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS/test2:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC1_rsfMRI/CSUB-{subject}C-01/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00014 00019 -ss 01 -c dcm2niix -b --overwrite

sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS/test2:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC2_T1/CSUB-{subject}C-02/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00014 00019 -ss 02 -c dcm2niix -b --overwrite

sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS/test2:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC2_rsfMRI/CSUB-{subject}C-02/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00014 00019 -ss 02 -c dcm2niix -b --overwrite


##
sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS/test2:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC1_*/CSUB-{subject}C-01/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00014 -ss 01 -c dcm2niix -b --overwrite ; sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS/test2:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC2_*/CSUB-{subject}C-02/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00014 -ss 02 -c dcm2niix -b --overwrite

##
sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS/test2:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC1_*/CSUB-{subject}C-01/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00014 00019 00026 -ss 01 -c dcm2niix -b --overwrite ; sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS/test2:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC2_*/CSUB-{subject}C-02/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00014 00019 00026 -ss 02 -c dcm2niix -b --overwrite


##
sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS/test3:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC1_*/CSUB-{subject}C-01/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00014 -ss 01 -c dcm2niix -b --overwrite ; sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS/test3:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC2_*/CSUB-{subject}C-02/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00014 -ss 02 -c dcm2niix -b --overwrite

##
sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS/test4:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC1_*/CSUB-{subject}C-01/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00014 00019 00026 -ss 01 -c dcm2niix -b --overwrite ; sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS/test4:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC2_*/CSUB-{subject}C-02/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00014 00019 00026 -ss 02 -c dcm2niix -b --overwrite


##
sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS/test_2sub:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC1_*/CSUB-{subject}C-01/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00014 00019 -ss 01 -c dcm2niix -b --overwrite

sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS/test_2sub:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC2_*/CSUB-{subject}C-02/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00014 00019 -ss 02 -c dcm2niix -b --overwrite



##
sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC1_*/CSUB-{subject}C-01/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00001 00004 00006 00007 00008 00009 00010 00011 00014 00015 00019 00020 00021 00023 00026 00028 00029 00030 00032 00035 00037 00040 00043 00044 00045 00046 00047 00048 00049 00050 00051 00052 00053 00055 00057 00060 00061 00062 00064 00065 00066 00067 00068 00069 00070 00071 00072 00074 00076 00079 00081 00084 00085 00086 00087 00088 00089 00090 00091 00093 00095 00097 00098 00099 00100 00101 00102 00104 00105 00107 00108 00110 00111 00112 00113 00115 00116 00119 00123 00124 00129 00130 00131 00132 00133 00136 00137 00138 00139 00140 00142 00143 00144 00145 00147 00148 00149 00150 00151 00152 00154 00155 00156 00157 00158 00160 00161 00162 00163 00164 00165 00166 00167 00168 00170 00171 00173 00174 00176 00177 00179 00180 00185 00186 00187 00188 00189 00190 00191 00192 00193 00195 00196 00198 00199 00200 00201 00202 00203 00204 00205 00206 00208 00210 00212 00213 00216 00220 00221 00222 00223 00225 00226 00228 00231 00232 00235 00236 00237 00238 00239 00241 00242 00244 00245 00246 00247 00248 00249 00250 00254 00255 00256 00257 00259 00262 00263 00264 00265 00266 00267 00268 00269 00270 00273 00275 00276 00277 00278 00279 00281 00282 00283 00284 00286 00288 00290 00291 00292 00293 00295 00296 00297 00298 00300 00301 00303 00304 00305 -ss 01 -c dcm2niix -b --overwrite ; sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC2_*/CSUB-{subject}C-02/*.dcm -o /base/Nifti/ -f /base/Nifti/code/heuristic.py -s 00003 00004 00006 00007 00009 00010 00011 00012 00014 00015 00018 00019 00020 00021 00022 00023 00024 00026 00027 00029 00031 00032 00037 00039 00040 00041 00043 00044 00045 00046 00047 00048 00049 00050 00051 00052 00053 00055 00056 00057 00060 00062 00065 00066 00068 00070 00071 00074 00076 00079 00084 00086 00087 00088 00089 00090 00091 00093 00095 00096 00097 00098 00099 00100 00101 00104 00105 00107 00109 00110 00111 00113 00115 00120 00123 00124 00129 00131 00132 00136 00137 00138 00142 00143 00144 00145 00147 00148 00149 00152 00154 00155 00156 00157 00158 00160 00161 00162 00163 00164 00165 00166 00167 00168 00171 00173 00174 00176 00177 00180 00183 00186 00187 00188 00189 00190 00191 00192 00194 00195 00196 00198 00199 00200 00201 00202 00203 00204 00205 00206 00207 00208 00210 00211 00212 00213 00216 00220 00221 00222 00223 00225 00226 00228 00230 00232 00237 00238 00239 00240 00241 00244 00245 00246 00248 00249 00250 00253 00254 00255 00256 00257 00258 00260 00263 00264 00265 00266 00269 00273 00275 00276 00277 00281 00284 00285 00288 00292 00293 00299 00300 00301 00302 00303 00306 00307 00308 00309 00310 00312 00313 00314 00315 00316 00318 00322 00323 00324 00325 00326 00327 00328 00329 00330 00332 00333 00334 00335 00336 00337 00338 00339 00341 00343 00344 00345 00346 00347 00348 00349 00350 00351 -ss 02 -c dcm2niix -b --overwrite