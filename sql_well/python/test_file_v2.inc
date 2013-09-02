# well declaration W_SPEC
12.04.2011 ; *; W_SPEC; p001; 100.4; 89
* ; 12:32:48.22; W_SPEC; i002;
*; *; W_SPEC; p003

# well branch W_BRANCH_F
*;*;W_BRANCH_F;p001;*;*;*;p001.dev
*;*;W_BRANCH_F;p003;*;*;*;p003.dev
*;*;W_BRANCH_F;i002;*;*;*;i002.dev

# well completions W_COMP
# date      ; time      ; id        ; well      ; branch    ; STATUS    ; md        ; l     ; r     ; skin  ; khmult
#---------------------------------------------------------------------------------------------------------------------
01.02.2011  ;*          ;W_COMP     ;p001       ;main       ;OPEN       ;2600       ; 100;
01.03.2011  ;*          ;W_COMP     ;p003       ;main       ;OPEN       ;2500       ;200    ;0.1    ;1.0    ;2.0   
01.04.2011  ;*          ;W_COMP     ;i002       ;*          ;*          ;2500       ;200    ;0.1    ;1.0   

# well frac W_FRAC
# date      ; time      ; id        ; well      ; branch    ; STATUS    ; md        ; angle     ; H_L_1     ; H_L_2     ; H_U   ; H_D   ; PERM  ; H_TH  ;
#--------------------------------------------------------------------------------------------------------------------------------------------------------
*           ;*          ;W_FRAC     ;p003       ;main       ;OPEN       ;2495       ; 46        ; 240       ; 200       ; 3     ; 5     ;
*           ;*          ;W_FRAC     ;p003       ;*          ;OPEN       ;2795       ; 46        ; 240       ; 200       ; 3     ; 5     ;
*           ;*          ;W_FRAC     ;p003       ;main       ;SHUT       ;2495       ; 46        ; 240       ; 200       ; 3     ; 5     ;



# well PROD
*;*;W_PROD;p001;*;*;50;
12.05.2010;*;W_PROD;p001;*;*;50;
*;*;W_PROD;p003;*;*;40;

#well INJ
*;*;W_INJ;i002;*;*;*;300;



