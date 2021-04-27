; ---------------------------------------------------------
; Regressione con istruzioni SSE a 32 bit
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
;     nasm -f elf32 regression32.nasm 
;
%include "sseutils32.nasm"

section .data			; Sezione contenente dati inizializzati
	align	16
	z	dd	0,0,0,0
	align	16
	eps	dd	0.00000001,0.00000001,0.00000001,0.00000001
section .bss			; Sezione contenente dati non inizializzati
	alignb 16
	w resd 1
	;alignb 16
	;v1 resd 4
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
	mov	eax, %1
	push	eax
	mov	eax, %2
	push	eax
	call	get_block
	add	esp, 8
%endmacro

%macro	fremem	1
	push	%1
	call	free_block
	add	esp, 4
%endmacro

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------

N	equ	16
B	equ	12
A	equ	8

SOMMA	equ 		8
X		equ		12
FATTORE	equ		16
DIM		equ		20
G		equ		24
SCALARE	equ		28


global prodottoScalare
global aggiornaTheta
global aggiornaThetaAdagrad

aggiornaThetaAdagrad:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
	
	push		ebp							
	mov			ebp, esp					
	push		ebx							
	push		esi
	push		edi
	
	; ------------------------------------------------------------
	; legge i parametri dal Record di Attivazione corrente
	; ------------------------------------------------------------
	MOV		EAX,[EBP+SOMMA]
	MOV		EBX,[EBP+X]
	MOV 		ECX,[EBP+G]
	MOV		ESI,[EBP+DIM]
	
	XOR			EDI,EDI
	XOR			EDX,EDX
	SHR			ESI,2
	
	MOVSS		XMM2,[EBP+FATTORE]		; eta/v
	SHUFPS		XMM2,XMM2,0x00
	MOVSS		XMM3,[EBP+SCALARE]		; prodottoScalare
	SHUFPS		XMM3,XMM3,0x00
	MOVAPS		XMM6,[eps]
	
c7:	MOVAPS		XMM0,[EAX+EDX]			; carico 4 float da theta
	;MOVAPS [v1],XMM0
        ;printps v1,1
	MOVUPS		XMM1,[EBX+EDX]			; carico 4 float da xast
	;MOVAPS [v1],XMM1
        ;printps v1,1
	MOVUPS		XMM5,[ECX+EDX]			; carico 4 float da G
	;MOVAPS [v1],XMM5
        ;printps v1,1
	MULPS		XMM1,XMM3				; pScalare * xi
	MOVAPS		XMM4,XMM1				; carico pScalare * xi
	MULPS		XMM1,XMM2				; eta/v * pScalare * xi
	
	MULPS		XMM4,XMM4				; gj^2
	ADDPS		XMM5,XMM4				; eseguo sommatoria Gj
	MOVUPS		[ECX+EDX],XMM5 			; salvo Gj
	
	ADDPS		XMM5,XMM6				;aggiungo eps a Gj
	SQRTPS		XMM5,XMM5				;eseguo radice della somma
			
	DIVPS		XMM1,XMM5				; eseguo divisione tra pScalre*x[i]*eta/v e radice di Gj+eps
	SUBPS		XMM0,XMM1				; aggiorno tetha
	
	MOVAPS		[EAX+EDX],XMM0			; lo salvo
	
	ADD EDX,16
	ADD EDI,1
	CMP EDI,ESI
	JL c7
	;ciclo per i restanti che esegue le stesse op
	SHL 		ESI,2
	MOV		ECX,[EBP+DIM]
	SUB			ECX,ESI
	XOR			EDI,EDI
	MOV		ESI,[EBP+G]
	JMP			if3
	
c8:	MOVSS		XMM0,[EAX+EDX]
	MOVSS		XMM1,[EBX+EDX]
	MOVSS		XMM5,[ESI+EDX]
	
	MULSS		XMM1,XMM3
	MOVSS		XMM4,XMM1
	MULSS		XMM1,XMM2
	
	MULSS		XMM4,XMM4
	ADDSS		XMM5,XMM4
	MOVSS		[ESI+EDX],XMM5
	ADDSS		XMM5,XMM6
	
	SQRTSS		XMM5,XMM5
	DIVSS		XMM1,XMM5
	SUBSS		XMM0,XMM1
	
	MOVSS		[EAX+EDX],XMM0
	
	ADD			EDX,4
	ADD			EDI,1
if3:	CMP		EDI,ECX
	JL			c8
	
	pop	edi									; ripristina i registri da preservare
	pop	esi
	pop	ebx
	
	mov	esp, ebp							; ripristina lo Stack Pointer
	
	pop	ebp
					; ripristina il Base Pointer
	ret	

aggiornaTheta:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
	
	push		ebp							
	mov			ebp, esp					
	push		ebx							
	push		esi
	push		edi
	
	; ------------------------------------------------------------
	; legge i parametri dal Record di Attivazione corrente
	; ------------------------------------------------------------
	
	MOV EAX, [EBP+SOMMA]
	MOV EBX, [EBP+X]
	MOV ESI,	[EBP+DIM]
	
	XOR EDI,EDI
	XOR EDX,EDX
	SHR ESI,2
	MOVSS XMM2,[EBP+FATTORE]
	SHUFPS XMM2,XMM2, 0x00		;duplico valore su tutte le 4 celle del registro
	;MOVAPS [v1],XMM2
        ;printps v1,1
c3:	MOVAPS	XMM0,[EAX+EDX] 		; carico 4 float dalla somma
	;MOVAPS [v1],XMM0
        ;printps v1,1
	MOVUPS	XMM1,[EBX+EDX]		; carico 4 float da xast
	;MOVAPS [v1],XMM1
        ;printps v1,1
	MULPS	XMM1,XMM2			; moltiplico xast per fattore
	;MOVAPS [v1],XMM1
        ;printps v1,1
	SUBPS	XMM0,XMM1			;aggiorno somma
	
	MOVAPS [EAX+EDX],XMM0		;aggiorno memoria
	ADD EDX,16
	ADD EDI,	1
	CMP EDI,ESI
	JL c3
	
	SHL ESI,2
	MOV ECX,[EBP+DIM]
	SUB ECX,ESI
	XOR EDI,EDI
	JMP	if1

c4:	MOVSS XMM0,[EAX+EDX]	
	MOVSS XMM1,[EBX+EDX]
	MULSS XMM1,[EBP+FATTORE]
	SUBSS XMM0,XMM1
	
	MOVSS [EAX+EDX],XMM0
	ADD EDX,4
	ADD EDI,1
if1:	CMP EDI,ECX
	JL c4
	
	pop	edi									; ripristina i registri da preservare
	pop	esi
	pop	ebx
	
	mov	esp, ebp							; ripristina lo Stack Pointer
	
	pop	ebp
					; ripristina il Base Pointer
	ret										; torna alla funzione C chiamante

prodottoScalare:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
	push		ebp							; salva il Base Pointer
	mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
	push		ebx							; salva i registri da preservare
	push		esi
	push		edi
	; ------------------------------------------------------------
	; legge i parametri dal Record di Attivazione corrente
	; ------------------------------------------------------------
	
	
	;MOVAPS XMM0,[z]
	xorps XMM0,XMM0
	
	XOR EDI,EDI
	XOR EDX,EDX
	MOV EAX, [EBP+A]
	MOV EBX,[EBP+B]
	MOV ESI,[EBP+N]
	SHR ESI,2
 c:    MOVAPS XMM1,[EAX+EDX] ;CARICO 4 FLOAT DA A
        
	;MOVAPS [v1],XMM1
        ;printps v1,1
        
        ;MOV EAX,[EBP+B]
        MOVUPS XMM2,[EBX+EDX]
        ;MOVUPS [v2],XMM2
        ;printps v2,1
        
        ;MOV EAX,[EBP+N]
        ;MOV [w],EAX
        ;printss w
        
        MULPS XMM1,XMM2 ;MOLTIPLICO FRA DI LORO I MEMBRI
        ADDPS XMM0,XMM1
        ADD EDX,16
	ADD EDI,1
	CMP EDI,ESI
	Jl c
	
	SHL ESI,2
	MOV ECX,[EBP+N]
	SUB ECX,ESI
	XOR EDI,EDI
JMP	if
c2:	MOVSS	XMM1,[EAX+EDX]
	MULSS	XMM1,[EBX+EDX]
	ADDPS 	XMM0,XMM1
	ADD EDX,4
	ADD EDI,1
if:	CMP EDI,ECX
	JL c2
	
	HADDPS XMM0,XMM0
        HADDPS XMM0,XMM0
        MOVSS [w], XMM0
        
	;printss w
	fld dword [w]
	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------

		
	pop	edi									; ripristina i registri da preservare
	pop	esi
	pop	ebx
	
	mov	esp, ebp							; ripristina lo Stack Pointer
	
	pop	ebp
					; ripristina il Base Pointer
	ret										; torna alla funzione C chiamante
