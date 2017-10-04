#include "license.h"

#define INTERFACE "eth0"

int ip_address(char is_me, char *address)
{
     struct ifreq ifr;
     struct sockaddr_in *sa = (struct sockaddr_in *)&ifr.ifr_addr;
     int sockfd;

     memset(&ifr, 0, sizeof ifr);
     if (0 > (sockfd = socket(AF_INET, SOCK_STREAM, 0))) { perror("socket()"); return EXIT_FAILURE; }
     strcpy(ifr.ifr_name, INTERFACE);
     sa->sin_family = AF_INET;
     if (0 > ioctl(sockfd, SIOCGIFADDR, &ifr)) { perror("ioctl()"); return EXIT_FAILURE; }
     // if(is_me) fprintf(stderr,"%s : [%s]\n", ifr.ifr_name, inet_ntoa(sa->sin_addr));
     sprintf(address,"%s",inet_ntoa(sa->sin_addr));
     return EXIT_SUCCESS;
}

Int4 mac_addr_sys(u_char *addr)
/* implementation for Linux */
{
    struct ifreq ifr;
    struct ifreq *IFR;
    struct ifconf ifc;
    char buf[1024];
    int s, i;
    int ok = 0;

    s = socket(AF_INET, SOCK_DGRAM, 0);
    if (s==-1) { return -1; }
    ifc.ifc_len = sizeof(buf);
    ifc.ifc_buf = buf;
    ioctl(s, SIOCGIFCONF, &ifc);
    IFR = ifc.ifc_req;
    for (i = ifc.ifc_len / sizeof(struct ifreq); --i >= 0; IFR++) {
        strcpy(ifr.ifr_name, IFR->ifr_name);
        if (ioctl(s, SIOCGIFFLAGS, &ifr) == 0) {
            if (! (ifr.ifr_flags & IFF_LOOPBACK)) {
                if (ioctl(s, SIOCGIFHWADDR, &ifr) == 0) { ok = 1; break; }
            }
        }
    }
    close(s);
    if(ok){ bcopy( ifr.ifr_hwaddr.sa_data, addr, 6); }
    else { return -1; }
    return 0;
}

#define LICENSE_STARTUP "Copyright 1997-2017 Cold Spring Harbor Laboratory\n\
and the University of Maryland School of Medicine.\n\n\
All rights reserved. In accordance with United States Code Title 17,\n\
no portion of this software may be reproduced in any form or by any \n\
means, without permission in writing from the Office of Technology \n\
Transfer at Cold Spring Harbor Laboratory and the University of Maryland School of Medicine.\n\n\
This software is distributed by Cold Spring Harbor Laboratory \n\
and the University of Maryland School of Medicine only to Licensed Users.\n\
For further information contact:\n\
 Andrew F. Neuwald\n\
 Institute for Genome Sciences and\n\
 Department of Biochemistry & Molecular Biology\n\
 University of Maryland School of Medicine\n\
 Baltimore, MD 21201\n\
 Tel: 410-706-6724; E-mail: aneuwald@som.umaryland.edu\n\
 Home page: http://www.medschool.umaryland.edu/profiles/Neuwald-Andrew\n\n"

#define LICENSE_STARTUP2 "Copyright 1997-2017 the University of Maryland School of Medicine.\n\n\
All rights reserved. In accordance with United States Code Title 17,\n\
no portion of this software may be reproduced in any form or by any \n\
means, without permission in writing from the Office of Technology \n\
Transfer at the University of Maryland School of Medicine.\n\n\
This software is distributed by the University of Maryland School of Medicine only to Licensed Users.\n\
For further information contact:\n\
 Andrew F. Neuwald\n\
 Institute for Genome Sciences and\n\
 Department of Biochemistry & Molecular Biology\n\
 University of Maryland School of Medicine\n\
 Baltimore, MD 21201\n\
 Tel: 410-706-6724; E-mail: aneuwald@som.umaryland.edu\n\
 Home page: http://www.medschool.umaryland.edu/profiles/Neuwald-Andrew\n\n"

#define LICENSE_STARTUP3 "    Copyright (C) 1997-2014 Andrew F. Neuwald, Cold Spring Harbor Laboratory\n\
    and the University of Maryland School of Medicine.\n\n\
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and \n\
    associated documentation files (the \"Software\"), to deal in the Software without restriction, \n\
    including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,\n\
    and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,\n\
    subject to the following conditions: \n\n\
    The above copyright notice and this permission notice shall be included in all copies or substantial\n\
    portions of the Software.\n\n\
    THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT\n\
    NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.\n\
    IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,\n\
    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE\n\
    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. \n\n\
    For further information contact:\n\
         Andrew F. Neuwald\n\
         Institute for Genome Sciences and\n\
         Department of Biochemistry & Molecular Biology\n\
         University of Maryland School of Medicine\n\
         Baltimore, MD 21201\n\
         Tel: 410-706-6724; Fax: 410-706-1482; E-mail: aneuwald@som.umaryland.edu\n\n"

#define LICENSE_STARTUP4 "  Copyright 2017 The University of Maryland\n\n\
 For information contact:\n\
   Andrew F. Neuwald\n\
   Institute for Genome Sciences\n\
   Department of Biochemistry & Molecular Biology\n\
   University of Maryland School of Medicine, Baltimore\n\
   E-mail: aneuwald@som.umaryland.edu\n"

static char LicenseWasPrinted='F';

void	TurnOffLicenseStatement(){ LicenseWasPrinted='T'; }

void	PrintLicenseStatement(char *program)
{	
	if(LicenseWasPrinted == 'F'){
	  fprintf(stderr,"===========================================================\n");
	  fprintf(stderr,"  %s %s",program,LICENSE_STARTUP4); LicenseWasPrinted='T'; 
	  fprintf(stderr,"===========================================================\n");
	}
}

void	PrintLicenseStatement()
// { if(LicenseWasPrinted == 'F') fprintf(stderr,LICENSE_STARTUP3); LicenseWasPrinted='T'; }
{ if(LicenseWasPrinted == 'F') fprintf(stderr,LICENSE_STARTUP2); LicenseWasPrinted='T'; }

char	IsLicenseExpired()
{
    Int4 stat;
    int i;
    char okay=0,is_me=0;
    char str[1002];
    static Int4 no_calls=0;

#if 1
    PrintLicenseStatement(); return 0;
#endif
    // fprintf(stderr,"call number %d\n",no_calls); 
  if(no_calls==0){
    no_calls++;
#if 0
    //********************* Check user name *******************
    okay=0;
    char *user=getenv("USER");
    // if(user && strcmp("aneuwald",user)!=0){
    if(user){
       if(strcmp("afn_guest1",user) == 0){ return 0; // skip for afn_guest1...
       } else if(strcmp("afn_guest2",user) == 0){ return 0; // skip for afn_guest2...
       } else PrintLicenseStatement();
    } // else return 0; // skip this for me...
    if(user){
	if(strcmp("aneuwald",user)==0){ okay=1; is_me=1; }
	if(strcmp("nkannan",user)==0) okay=1;
	if(strcmp("kannan",user)==0) okay=1;
	if(strcmp("amar",user)==0) okay=1;
	if(strcmp("esbg",user)==0) okay=1;
	if(strcmp("etal",user)==0) okay=1;
	// if(is_me) fprintf(stderr,"user: %s\n",user);
    } else print_error("fatal: licensing error 1.");
    if(!okay) print_error("Fatal: unlicensed user.");

    //********************* Check host name *******************
    okay=0;
    char host[256];
    int status=gethostname(host,255);
    if(status==0){
        if(is_me) fprintf(stderr,"host: %s\n",host);
	if(strcmp("bullseye.igs.umaryland.edu",host)==0) okay=1;
	if(strcmp("aneuwald-lx.igs.umaryland.edu",host)==0) okay=1;

	// Kannan at Georgia:
	if(strcmp("rcluster.rcc.uga.edu",host)==0) okay=1;
	if(strcmp("esbg-work1",host)==0) okay=1;
	if(strcmp("esbg-work1.rcc.uga.edu",host)==0) okay=1;
	if(strcmp("esbg-work2",host)==0) okay=1;
	if(strcmp("esbg-work2.rcc.uga.edu",host)==0) okay=1;
	if(strcmp("esbg-work3",host)==0) okay=1;
	if(strcmp("esbg-work3.rcc.uga.edu",host)==0) okay=1;

	if(strcmp("kannan-wireless",host)==0) okay=1;
	if(strcmp("KANNAN",host)==0) okay=1;
    } else print_error("fatal: licensing error 2.");
    if(!okay) print_error("Fatal: unlicensed machine (1).");
#endif

#if 0
    //********************* Check IP address *******************
    okay=0;		// /sbin/ifconfig
    status=ip_address(is_me,str);
    if(status==0){
	if(strcmp("143.48.42.3",str)==0) okay=1; // photon
	if(strcmp("143.48.42.5",str)==0) okay=1; // electron
	if(strcmp("143.48.1.75",str)==0) okay=1; // quark

	if(strcmp("134.192.146.172",str)==0) okay=1; // aneuwald-lx.igs.umaryland.edu
	// Kannan at Georgia:
	if(strcmp("128.192.75.100",str)==0) okay=1; // quark

	if(strcmp("143.48.42.128",str)==0) okay=1; // baryon
	if(strcmp("132.239.16.141",str)==0) okay=1; // livingstone.
	if(strcmp("132.239.164.85",str)==0) okay=1; // tioga
	if(strcmp("132.239.255.255",str)==0) okay=1; // tioga
	if(strcmp("143.48.7.47",str)==0) okay=1; // kannan-wireless
	if(strcmp("132.239.164.94",str)==0) okay=1; // sanmarcos
    } else print_error("fatal: licensing error 3.");
    if(!okay) print_error("Fatal: unlicensed machine (2).");

    //********************* Check CPU properties *******************
#endif
#if 0
    okay=0;
    FILE *fp=fopen("/proc/cpuinfo","r");
    if(fp==NULL) print_error("Fatal: licensing error 3");
    while(fgets(str,1000,fp) != NULL){
	// note: need two tab after 'MHz'.
	// model name      : Intel(R) Xeon(TM) CPU 3.06GHz
	if(strcmp("cpu MHz		: 3065.902\n",str) ==0 ||
	   strcmp("model name	: Intel(R) Xeon(TM) CPU 3.06GHz\n",str) ==0) // electron.
		{ okay=1; break; }

	// Kannan at Georgia:
	if(strcmp("cpu MHz		: 3065.897\n",str) ==0 ||
	   strcmp("model name	: Dual Core AMD Opteron(tm) Processor 275\n",str) ==0) // rcluster
		{ okay=1; break; }

	if(strcmp("cpu MHz		: 3065.897\n",str) ==0 ||
	   strcmp("model name	: Intel(R) Xeon(TM) CPU 3.06GHz\n",str) ==0) // photon.
		{ okay=1; break; }
	if(strcmp("cpu MHz		: 2199.791\n",str) ==0 ||
	   strcmp("model name	: Intel(R) XEON(TM) CPU 2.20GHz\n",str) ==0) // quark
		{ okay=1; break; }
	if(strcmp("cpu MHz		: 2191.204\n",str) ==0 ||
	   strcmp("model name	: Intel(R) Xeon(TM) CPU 2.20GHz\n",str) ==0) // baryon
		{ okay=1; break; }
	if(strcmp("cpu MHz		: 931.316\n",str) ==0 ||
	   strcmp("model name	: Pentium III (Coppermine)\n",str) ==0) // livingstone.
		{ okay=1; break; }
	if(strcmp("cpu MHz		: 3200.160\n",str) ==0 ||
	   strcmp("model name	: Intel(R) Xeon(TM) CPU 3.20GHz\n",str) ==0) // tioga
		{ okay=1; break; }
	if(strcmp("cpu MHz		: 599.599\n",str) ==0 ||
	   strcmp("model name	: Intel(R) Pentium(R) M processor 1.80GHz\n",str)==0) // kannan
		{ okay=1; break; }
	if(strcmp("cpu MHz		: 3400.267\n",str) ==0 ||
	   strcmp("model name	: Intel(R) Xeon(TM) CPU 3.40GHz\n",str)==0) // sanmarcos
		{ okay=1; break; }
    } fclose(fp);
    if(!okay) print_error("Fatal: unlicensed machine (3).");
    // else if(is_me) fprintf(stderr,"%s",str);
#endif

#if 0
    //********************* Check MAC address *******************
    okay=0;
    u_char addr[6];
    stat = mac_addr_sys(addr);
    if(0 == stat) {
	sprintf(str,"%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x",
		addr[0],addr[1],addr[2],addr[3],addr[4],addr[5]);
        // if(is_me) fprintf(stderr,"MAC address: %s\n",str);
	if(strcmp("00096bb5edc6",str)==0) okay=1; // photon.cshl.edu
	if(strcmp("000d604e2702",str)==0) okay=1; // electron.cshl.edu
	if(strcmp("000476f05dae",str)==0) okay=1; // quark.cshl.edu
	if(strcmp("000e0c314628",str)==0) okay=1; // baryon.cshl.edu
	if(strcmp("00b0d03dcb32",str)==0) okay=1; // livingstone
	if(strcmp("0011093b3ce5",str)==0) okay=1; // tioga
	if(strcmp("000b7d135b8d",str)==0) okay=1; // kannan
	if(strcmp("00145e45044e",str)==0) okay=1; // sanmarcos

	if(strcmp("001aa0ad75c0",str)==0) okay=1; // aneuwald-lx.igs.umaryland.edu
	// Kannan at Georgia
	if(strcmp("00e0814281d4",str)==0) okay=1; // kannan
    } else {
        print_error("Fatal: licensing error 4.");
        // fprintf(stderr,"can't get MAC address (status = %d)\n",stat); exit( 1);
    }
    if(!okay) print_error("Fatal: unlicensed machine.");
    // else fprintf(stderr,"MAC address: %s\n",str);
#endif

    //********************* Check to see whether license has expired *******************
    Int4 start_time = 1455891521; // Fri Feb 19 09:19:17 EST 2016
#if 0
    Int4 start_time = 1441726180; // Tue Sep  8 11:30:24 EDT 2015
    Int4 start_time = 1437412974; // Mon Jul 20 13:23:09 EDT 2015
    Int4 start_time = 1431975002; // Mon May 18 14:50:28 EDT 2015
    Int4 start_time = 1427469366; // Fri Mar 27 11:17:03 EDT 2015
    Int4 start_time = 1420643066; // Wed Jan  7 10:05:22 EST 2015
    // March 31 2005

    // about May 28, 2006...= (+2 years == 63072000)
    Int4 expiration_time = start_time+63072000;
    // Feb. 3, 2005
    Int4 start_time = 1107439911;
#else

    //  2 years == 63072000
    // 400 days == 34560000)
    Int4 expiration_time = start_time+34560000; // add 400 days
    // Int4 expiration_time = start_time+(34560000/2); // add 200 days
    // Int4 expiration_time = start_time+7776000; // add 90 days
    // Int4 expiration_time = start_time+3456000; // add 40 days
    // Int4 expiration_time = start_time; // test
#endif

    Int4 current_time=time(NULL);
    // fprintf(stderr,"%d (%d)\n",current_time,expiration_time);
    // fprintf(stderr,"%d days used/",(current_time-start_time)/86400);
    // if(strcmp("aneuwald",user)==0){ okay=1; is_me=1; }
    // PrintLicenseStatement();
    fprintf(stderr,"%d days left on license.\n",(expiration_time-current_time)/86400);

    if(current_time > expiration_time){
           print_error("License expired: contact aneuwald@som.umaryland.edu.");
    }
  } return 0;
}

