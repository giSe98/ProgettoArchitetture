; ---------------------------------------------------------
; Regression con istruzioni AVX a 64 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mediante il packaging tool del sistema 
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;     nasm -f elf64 regression64.nasm
;

%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati
align	32
eps	dq	0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001
section .bss			; Sezione contenente dati non inizializzati

alignb 32
eta		resq	1
alignb 32
w		resq 4
section .text			; Sezione contenente il codice macchina

; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block

%macro	getmem	2
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro

global aggiornaTheta
global aggiornaThetaAdagrad
global prodottoScalare

;msg	db 'eta:',0
;nl	db 10,0

prodottoScalare:
		;sequenza ingresso
        push		rbp				; salva il Base Pointer
        mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
        pushaq		; salva i registri general

	MOV 	RAX,RDI	;indirizzo A
	MOV	RBX,RSI	;indirizzo B
	MOV 	RSI,RDX	;valore di n
	SHR		RSI,2
	XOR 	RDI,RDI
	XOR		R8,R8
	MOV	R9,RCX ;indirizzo c
       VXORPD 	YMM0, YMM0 ;SOMMA=0
c1:   
        VMOVAPD  YMM1, [RAX+R8] ;copio 4 elementi di a
	;VMOVAPD	[w], YMM1
	;printpd	w,2

        ;VMOVUPD YMM2, [RBX+R8] ; copio 4 elementi di b
	;VMOVAPD	[w], YMM2
	;printpd	w,2

        VMULPD  	YMM1, [RBX+R8] ; a*b
	;VMOVAPD	[w], YMM1
	;printpd	w,2

        VADDPD  YMM0, YMM1 ; somma=somma+a*b
	;VMOVAPD	[w], YMM0
	;printpd	w,2

	ADD 	R8, 32 ;sposto il puntatore di a di 4 elementi
        ADD 	RDI, 1 
        CMP 	RDI,RSI 
        JL 		c1
	
	SHL 	RSI,2
	MOV	RCX,RDX
	SUB		RCX,RSI
	XOR		RDI,RDI
	JMP	if
c2: 	VMOVSD	XMM1,[RAX+R8]	; carico i residui di a
	;VMOVAPD	[w], YMM1
	;printpd	w,2
	;VMOVSD	XMM2,[RBX+R8]	; carico i residui di b
	;VMOVAPD	[w], YMM2
	;printpd	w,2
	
	VMULSD	XMM1,[RBX+R8]		; residui di a * residui di b
	;VMOVAPD	[w], YMM1
	;printpd	w,2
	VMOVSD XMM2,XMM1
	
	VADDPD	YMM0,YMM2		
	;VMOVAPD	[w], YMM0
	;printpd	w,2
	
	ADD		R8,8
	ADD		RDI,1
if:	CMP 	RDI,RCX
	JL		c2
	
	VHADDPD	YMM0,YMM0
	;VHADDPD	YMM0,YMM0
	;VMOVAPD	[w], YMM0
	;printpd	w,2
	
	VMOVAPD	YMM1, YMM0
	VPERM2F128	YMM1, YMM1, YMM1, 00000001b
	;VMOVAPD	[w], YMM1
	;printpd	w,2
	
	;VBLENDPD	YMM0, YMM1, 0100b
		
	;VMOVAPD	[w], YMM0
	;printpd	w,2
	VADDPD	YMM0,YMM1
	;VMOVAPD	[w], XMM0
	;printpd	w,1
	
	;VEXTRAXTPS	RAX,XMM0,0
	;VMOVAPD	[w], RAX
	;printpd	w,1
	
	VMOVSD	[R9],XMM0
	
	;sequenza uscita
        popaq						; ripristina i registri generali
        mov		rsp, rbp			; ripristina lo Stack Pointer
        pop		rbp					; ripristina il Base Pointer
        ret							; torna alla funzione C chiamante


aggiornaTheta:
        ;sequenza ingresso
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri general

		;fattore è un double quindi lo trovo in xmmo, non posso toccarlo
		;rdi theta allineato
		;rsi x non allineato
		;rdx dim
		VBROADCASTSD YMM2, XMM0 ; COPIO 4 VOLTE FATTORE IN YMM2
		;VMOVAPD	[w], YMM2
		;printpd	w,2
		XOR		R9,R9
		XOR		R8,R8
		MOV	RAX,RDX
		SHR 	RAX,2


c6:		VMOVAPD YMM0, [RDI+R9]
		VMOVUPD YMM1, [RSI+R9]
		VMULPD 	YMM1, YMM2
		VSUBPD 	YMM0, YMM1
		VMOVAPD [RDI+R9], YMM0
		
		ADD R9, 32
		ADD R8,1
		CMP R8, RAX
		JL c6
		
		SHL	RAX,2
		SUB RDX,RAX						;	R9 = EDX		R8 = EDI		RAX	=  ESI	RDX = ECX
		XOR	R8,R8
		JMP if1
		
c7:		VMOVSD	XMM0, [RDI+R9]
		VMOVSD 	XMM1, [RSI+R9]
		VMULSD	XMM1,XMM2
		
		VSUBSD	XMM0,XMM1
		VMOVSD	[RDI+R9],XMM0
		ADD R9,8
		ADD R8,1
if1:		CMP R8,RDX
		JL	c7
		
		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp					; ripristina il Base Pointer
		ret							; torna alla funzione C chiamante


aggiornaThetaAdagrad:
	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

	;RDI =theta
	;RSI= x
	;RCX= g
	;RDX =dim
	;XMMO= fattore
	;XMM1= pscalare

	VBROADCASTSD YMM2, XMM0 ;YMM2 copio 4 volte fattore
	VBROADCASTSD YMM3, XMM1; YMM3 COPIO 4 VOLTE PS
	VMOVAPD YMM6, [eps]
	XOR R8,R8  ;CORRISPONDE A EDX
	XOR R9,R9  ;CORRISPONDE A EDI
	MOV RAX,RDX  ; 
	SHR RAX,2  ; CORRISPONDE A ESI
c8:
	VMOVAPD YMM0, [RDI+R8] ;THETA
	VMOVUPD YMM1, [RSI+R8] ;XAST
	VMOVUPD YMM5, [RCX+R8] ;GUARDARE SE G È ALLINEATO
	
	VMULPD YMM1, YMM3; PS*XI
	VMOVAPD YMM4, YMM1
	VMULPD YMM1, YMM2

	VMULPD YMM4, YMM4 ;gj^2
	VADDPD YMM5, YMM4 ;SOMMATORIA GJ
	VMOVUPS [RCX+R8], YMM5; SALVO GJ

	VADDPD YMM5, YMM6;EPS+GJ
	VSQRTPD YMM5, YMM5 ; RADICE DELLA SOMMA
	
	VDIVPD YMM1, YMM5 ; eseguo divisione tra pScalre*x[i]*eta/v e radice di Gj+eps
	VSUBPD YMM0, YMM1 ; AGGIORNO THETA
	
	VMOVAPS [RDI+R8], YMM0
	ADD R8, 32
	ADD R9, 1
	CMP R9,RAX
	JL c8
	;CICLO RESTANTI ELEMENTI
	SHL RAX,2
	SUB RDX,RAX
	XOR R9,R9
	JMP if3
	
c9:	VMOVSD XMM0, [RDI+R8]
	VMOVSD XMM1, [RSI+R8]
	VMOVSD XMM5, [RCX+R8]
	
	VMULSD XMM1, XMM3
	VMOVSD XMM4, XMM1
	VMULSD XMM1, XMM2
	
	VMULSD XMM4, XMM4
	VADDSD XMM5, XMM4
	VMOVSD [RCX+R8], XMM5
	VADDSD XMM5, XMM6
	
	VSQRTSD XMM5, XMM5
	VDIVSD XMM1, XMM5
	VSUBSD XMM0, XMM1
	
	VMOVSD [RDI+R8], XMM0
	
	ADD R8, 8
	ADD R9, 1
if3:	CMP R9,RDX
	JL	c9
	
	;sequenza uscita
        popaq						; ripristina i registri generali
        mov		rsp, rbp			; ripristina lo Stack Pointer
        pop		rbp					; ripristina il Base Pointer
        ret	
