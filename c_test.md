```c
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include <gsl/gsl_cdf.h>

#define SAMPLE_MARK "!Sample_title"
//271
//Group_1=53,Group_2=120,Group_3=98
//gene_num=11098
struct sample{
	char sample_id[10];
	int age;
};
struct gene_info{
	char gene_id[20];
	double Group_1[53];
	double Group_2[120];
	double Group_3[98];
	double Group_1_mean_value,Group_2_mean_value,Group_3_mean_value;
	double Group_1_var,Group_2_var,Group_3_var;
	double mean1_2,mean2_3,mean1_3;
	double var1_2,var2_3,var1_3;
};
struct t_test{
	char passed_gene_id[20];
	double p_value;
};
void swap_t_test(struct t_test *a,struct t_test *b){
	struct t_test temp = *a;
	*a = *b;
	*b = temp;
}
void sort_t_test_by_p(struct t_test *arr,int n){
	for(int i=0;i<n-1;i++){
		for(int j=0;j<n-i-1;j++){
			if(arr[j].p_value>arr[j+1].p_value)
				swap_t_test(&arr[j],&arr[j+1]);
		}
	}
}
int min_get(int cnt){
	if(cnt<50)
		return cnt;
	if(cnt>=50){
		return 50;
	}
}
//count gene numbers 
int gene_count(){
	int gene_num=0;
	char line[200];
	FILE *p1 = fopen("GSE277909_series_matrix.txt","r");
	while(fgets(line,sizeof(line),p1)!=NULL){
		if(line[0]=='E'&&line[1]=='N'&&line[2]=='S'&&line[3]=='G')
			gene_num++;
	}
	fclose(p1);
	return gene_num;
}
//count sample numbers 
int lines_count(){
	FILE *p1 = fopen("GSE277909_ID_Age.txt","r");
	if(p1==NULL){
		printf("Can Not open the file");
		return -1;
	}
	char line[100];
	int line_count=0;
	if(fgets(line,sizeof(line),p1)==NULL){
		fclose(p1);
		return 0;
	}
	while(fgets(line,sizeof(line),p1)!=NULL)
		line_count++;
	fclose(p1);
	return line_count;
}
//get age into struct sample
struct sample * id_age_get(){
	int line_num=lines_count();
	char line[100];
	struct sample *sample_arr=(struct sample *)malloc((line_num*(sizeof(struct sample))));
	FILE *p1 = fopen("GSE277909_ID_Age.txt","r");
	if(p1==NULL){
		printf("Can Not open the file");
		return NULL;
	}
	if(fgets(line,sizeof(line),p1)==NULL){
		fclose(p1);
		return 0;
	}
	for(int i=0;i<line_num;i++){
		fscanf(p1," %s",sample_arr[i].sample_id);
		fscanf(p1," %d",&sample_arr[i].age);
	}
	return sample_arr;
}
int * id_order_in_matrix(struct sample *p){
	FILE *p1 = fopen("GSE277909_series_matrix.txt","r");
	char no_use_string[100];
	char id[271][10];
	while(1){
		fscanf(p1," %s",no_use_string);
		if(strcmp(no_use_string,SAMPLE_MARK)==0)
			break;
	}
	for(int i=0;i<271;i++){
		fscanf(p1," \"%[^\"]\" ",id[i]);
	}
	fclose(p1);
	int * order_arr=(int *)malloc(271*sizeof(int));
	for(int i=0;i<271;i++){
		for(int j=0;j<271;j++){
			if(strcmp((p+j)->sample_id,id[i])==0){
				if((p+j)->age<30){
					order_arr[i] = 1;
				}
				else if((p+j)->age>=30&&(p+j)->age<=50){
					order_arr[i] = 2;
				}
				else{
					order_arr[i] = 3;
				}
			}
		}
	}
	return order_arr;
}
struct gene_info * gene_info_get(int * group_order,int gene_nums){
	FILE *p1 = fopen("GSE277909_series_matrix.txt","r");
	struct gene_info * gene_arr = (struct gene_info *)malloc(gene_nums*sizeof(struct gene_info));
	char no_use_string[100];
	double sum_1=0,sum_2=0,sum_3=0;
	while(1){
		fscanf(p1," %s",no_use_string);
		if(no_use_string[0]=='E'&&no_use_string[1]=='N'&&no_use_string[2]=='S'&&no_use_string[3]=='G')
			break;
	}
	strcpy(gene_arr->gene_id,no_use_string);
	int a=0,b=0,c=0;
	for(int i=0;i<271;i++){
		if(*(group_order+i)==1){
			fscanf(p1," %lf",&(gene_arr->Group_1[a]));
			gene_arr->Group_1[a] = log2((gene_arr->Group_1[a])+1);
			sum_1+=gene_arr->Group_1[a];
			a++;
		}
		else if(*(group_order+i)==2){
			fscanf(p1," %lf",&(gene_arr->Group_2[b]));
			gene_arr->Group_2[b] = log2((gene_arr->Group_2[b])+1);
			sum_2+=gene_arr->Group_2[b];
			b++;
		}
		else{
			fscanf(p1," %lf",&(gene_arr->Group_3[c]));
			gene_arr->Group_3[c] = log2((gene_arr->Group_3[c])+1);
			sum_3+=gene_arr->Group_3[c];
			c++;
		}	
	}
	gene_arr->Group_1_mean_value = sum_1/a;	
	gene_arr->Group_2_mean_value = sum_2/b;
	gene_arr->Group_3_mean_value = sum_3/c;
	for(int i=1;i<gene_nums;i++){
		a=0,b=0,c=0;
		sum_1=0,sum_2=0,sum_3=0;
		fscanf(p1," %s",(gene_arr+i)->gene_id);	
		for(int j=0;j<271;j++){
			if(*(group_order+j)==1){
				fscanf(p1," %lf",&((gene_arr+i)->Group_1[a]));
				(gene_arr+i)->Group_1[a] = log2(((gene_arr+i)->Group_1[a])+1);
				sum_1+=(gene_arr+i)->Group_1[a];
				a++;
			}
			else if(*(group_order+j)==2){
				fscanf(p1," %lf",&((gene_arr+i)->Group_2[b]));
				(gene_arr+i)->Group_2[b] = log2(((gene_arr+i)->Group_2[b])+1);
				sum_2+=(gene_arr+i)->Group_2[b];
				b++;
			}
			else{
				fscanf(p1," %lf",&((gene_arr+i)->Group_3[c]));
				(gene_arr+i)->Group_3[c] = log2(((gene_arr+i)->Group_3[c])+1);
				sum_3+=(gene_arr+i)->Group_3[c];
				c++;
			}	
		}
		(gene_arr+i)->Group_1_mean_value = sum_1/a;	
		(gene_arr+i)->Group_2_mean_value = sum_2/b;
		(gene_arr+i)->Group_3_mean_value = sum_3/c;
	}
	return gene_arr;
}
void var_count(struct gene_info* gene_arr,int gene_nums){
	//double Group_1_var,Group_2_var,Group_3_var;
	double sum_1=0,sum_2=0,sum_3=0;
	for(int i=0;i<gene_nums;i++){
		sum_1=0,sum_2=0,sum_3=0;
		for(int a=0;a<53;a++){
			sum_1+=(((gene_arr+i)->Group_1[a])-((gene_arr+i)->Group_1_mean_value))*(((gene_arr+i)->Group_1[a])-((gene_arr+i)->Group_1_mean_value));
		}
		(gene_arr+i)->Group_1_var=sum_1/52;
		for(int b=0;b<120;b++){
			sum_2+=(((gene_arr+i)->Group_2[b])-((gene_arr+i)->Group_2_mean_value))*(((gene_arr+i)->Group_2[b])-((gene_arr+i)->Group_2_mean_value));
		}
		(gene_arr+i)->Group_2_var=sum_2/119;
		for(int c=0;c<98;c++){
			sum_3+=(((gene_arr+i)->Group_3[c])-((gene_arr+i)->Group_3_mean_value))*(((gene_arr+i)->Group_3[c])-((gene_arr+i)->Group_3_mean_value));
		}
		(gene_arr+i)->Group_3_var=sum_3/97;
	}
}
//count Group1+2 Group2+3 Group1+3
void mean_var_group_merge(struct gene_info* gene_arr,int gene_nums){
	double sum_1=0,sum_2=0,sum_3=0;
	for(int i=0;i<gene_nums;i++){
		(gene_arr+i)->mean1_2=(((gene_arr+i)->Group_1_mean_value)*53+((gene_arr+i)->Group_2_mean_value)*120)/(53+120);
		(gene_arr+i)->mean2_3=(((gene_arr+i)->Group_2_mean_value)*120+((gene_arr+i)->Group_3_mean_value)*98)/(120+98);
		(gene_arr+i)->mean1_3=(((gene_arr+i)->Group_1_mean_value)*53+((gene_arr+i)->Group_3_mean_value)*98)/(53+98);
	}
	for(int i=0;i<gene_nums;i++){
		sum_1=0,sum_2=0,sum_3=0;
		for(int a=0;a<53;a++){
			sum_1+=(((gene_arr+i)->Group_1[a])-((gene_arr+i)->mean1_2))*(((gene_arr+i)->Group_1[a])-((gene_arr+i)->mean1_2));
		}
		for(int b=0;b<120;b++){
			sum_1+=(((gene_arr+i)->Group_2[b])-((gene_arr+i)->mean1_2))*(((gene_arr+i)->Group_2[b])-((gene_arr+i)->mean1_2));
		}
		(gene_arr+i)->var1_2=sum_1/((53+120)-1);
		for(int b=0;b<120;b++){
			sum_2+=(((gene_arr+i)->Group_2[b])-((gene_arr+i)->mean2_3))*(((gene_arr+i)->Group_2[b])-((gene_arr+i)->mean2_3));
		}
		for(int c=0;c<98;c++){
			sum_2+=(((gene_arr+i)->Group_3[c])-((gene_arr+i)->mean2_3))*(((gene_arr+i)->Group_3[c])-((gene_arr+i)->mean2_3));
		}
		(gene_arr+i)->var2_3=sum_2/((120+98)-1);
		for(int a=0;a<53;a++){
			sum_3+=(((gene_arr+i)->Group_1[a])-((gene_arr+i)->mean1_3))*(((gene_arr+i)->Group_1[a])-((gene_arr+i)->mean1_3));
		}
		for(int c=0;c<98;c++){
			sum_3+=(((gene_arr+i)->Group_3[c])-((gene_arr+i)->mean1_3))*(((gene_arr+i)->Group_3[c])-((gene_arr+i)->mean1_3));
		}
		(gene_arr+i)->var1_3=sum_3/((53+98)-1);
	}
	
}
int main(){
	int line_num=lines_count(),Group_1_num=0,Group_2_num=0,Group_3_num=0,gene_nums;
	struct sample* sample_arr=id_age_get();
	struct sample *p;
	struct gene_info * gene_arr;
	p = sample_arr;
	gene_nums = gene_count();
	//count Group1 Group2 Group3 numbers to create struct gene_info
	for(int i=0;i<line_num;i++){
		if(((p+i)->age)<30)
			Group_1_num++;
	}
	for(int i=0;i<line_num;i++){
		if(((p+i)->age)>=30&&((p+i)->age)<=50)
			Group_2_num++;
	}
	for(int i=0;i<line_num;i++){
		if(((p+i)->age)>50)
			Group_3_num++;
	}
	//get the order of matrix
	int * order_arr = id_order_in_matrix(p);
	gene_arr = gene_info_get(order_arr,gene_nums);
	var_count(gene_arr,gene_nums);
	//Group1 VS Group2
	printf("Group1 VS Group2\n\n");
	struct t_test* t_test1=(struct t_test*)malloc(gene_nums*sizeof(struct t_test));
	double p1_test,term_1_1,term_1_2,t1,term_1_3,term_1_4,df1;
	int cnt1=0;
	for(int i=0;i<gene_nums;i++){
		term_1_1=(gene_arr+i)->Group_1_mean_value-(gene_arr+i)->Group_2_mean_value;
		term_1_2=pow((((gene_arr+i)->Group_1_var)/53)+(((gene_arr+i)->Group_2_var)/120),0.5);
		t1 = term_1_1/term_1_2;
		term_1_3=pow((((gene_arr+i)->Group_1_var)/53)+(((gene_arr+i)->Group_2_var)/120),2);
		term_1_4=(pow(((gene_arr+i)->Group_1_var)/53,2)/52)+(pow(((gene_arr+i)->Group_2_var)/120,2)/119);
		df1=term_1_3/term_1_4;
		p1_test=2*gsl_cdf_tdist_Q(fabs(t1),df1);
		if(p1_test<0.05){
			(t_test1+cnt1)->p_value=p1_test;
			strcpy((t_test1+cnt1)->passed_gene_id,(gene_arr+i)->gene_id);
			cnt1++;
		}
	}
	sort_t_test_by_p(t_test1,cnt1);
	int min1=min_get(cnt1);
	FILE *q1=fopen("Group1_VS_Group2.txt","w");
	for(int i=0;i<min1;i++){
		printf("GeneID:%s\tp-value:%.10lf\n",(t_test1+i)->passed_gene_id,(t_test1+i)->p_value);
		fprintf(q1,"GeneID:%s\tp-value:%.10lf\n",(t_test1+i)->passed_gene_id,(t_test1+i)->p_value);
	}
	fclose(q1);
	printf("\n\n");
	//Group2 VS Group3
	printf("Group2 VS Group3\n\n");
	struct t_test* t_test2=(struct t_test*)malloc(gene_nums*sizeof(struct t_test));
	double p2_test,term_2_1,term_2_2,t2,term_2_3,term_2_4,df2;
	int cnt2=0;
	for(int i=0;i<gene_nums;i++){
		term_2_1=(gene_arr+i)->Group_2_mean_value-(gene_arr+i)->Group_3_mean_value;
		term_2_2=pow((((gene_arr+i)->Group_2_var)/120)+(((gene_arr+i)->Group_3_var)/98),0.5);
		t2 = term_2_1/term_2_2;
		term_2_3=pow((((gene_arr+i)->Group_2_var)/120)+(((gene_arr+i)->Group_3_var)/98),2);
		term_2_4=(pow(((gene_arr+i)->Group_2_var)/120,2)/119)+(pow(((gene_arr+i)->Group_3_var)/98,2)/97);
		df2=term_2_3/term_2_4;
		p2_test=2*gsl_cdf_tdist_Q(fabs(t2),df2);
		if(p2_test<0.05){
			(t_test2+cnt2)->p_value=p2_test;
			strcpy((t_test2+cnt2)->passed_gene_id,(gene_arr+i)->gene_id);
			cnt2++;
		}
	}
	sort_t_test_by_p(t_test2,cnt2);
	int min2=min_get(cnt2);
	FILE *q2=fopen("Group2_VS_Group3.txt","w");
	for(int i=0;i<min2;i++){
		printf("GeneID:%s\tp-value:%.10lf\n",(t_test2+i)->passed_gene_id,(t_test2+i)->p_value);
		fprintf(q2,"GeneID:%s\tp-value:%.10lf\n",(t_test2+i)->passed_gene_id,(t_test2+i)->p_value);
	}
	fclose(q2);
	printf("\n\n");
	//Group1 VS Group3
	printf("Group1 VS Group3\n\n");
	struct t_test* t_test3=(struct t_test*)malloc(gene_nums*sizeof(struct t_test));
	double p3_test,term_3_1,term_3_2,t3,term_3_3,term_3_4,df3;
	int cnt3=0;
	for(int i=0;i<gene_nums;i++){
		term_3_1=(gene_arr+i)->Group_1_mean_value-(gene_arr+i)->Group_3_mean_value;
		term_3_2=pow((((gene_arr+i)->Group_1_var)/53)+(((gene_arr+i)->Group_3_var)/98),0.5);
		t3 = term_3_1/term_3_2;
		term_3_3=pow((((gene_arr+i)->Group_1_var)/53)+(((gene_arr+i)->Group_3_var)/98),2);
		term_3_4=(pow(((gene_arr+i)->Group_1_var)/53,2)/52)+(pow(((gene_arr+i)->Group_3_var)/98,2)/97);
		df3=term_3_3/term_3_4;
		p3_test=2*gsl_cdf_tdist_Q(fabs(t3),df3);
		if(p3_test<0.05){
			(t_test3+cnt3)->p_value=p3_test;
			strcpy((t_test3+cnt3)->passed_gene_id,(gene_arr+i)->gene_id);
			cnt3++;
		}
	}
	sort_t_test_by_p(t_test3,cnt3);
	int min3=min_get(cnt3);
	FILE *q3=fopen("Group1_VS_Group3.txt","w");
	for(int i=0;i<min3;i++){
		printf("GeneID:%s\tp-value:%.10lf\n",(t_test3+i)->passed_gene_id,(t_test3+i)->p_value);
		fprintf(q3,"GeneID:%s\tp-value:%.10lf\n",(t_test3+i)->passed_gene_id,(t_test3+i)->p_value);
	}
	fclose(q3);
	printf("\n\n");
	mean_var_group_merge(gene_arr,gene_nums);
	//Group1+2 VS Group3
	printf("Group1+2 VS Group3\n\n");
	struct t_test* t_test4=(struct t_test*)malloc(gene_nums*sizeof(struct t_test));
	double p4_test,term_4_1,term_4_2,t4,term_4_3,term_4_4,df4;
	int cnt4=0;
	for(int i=0;i<gene_nums;i++){
		term_4_1=(gene_arr+i)->mean1_2-(gene_arr+i)->Group_3_mean_value;
		term_4_2=pow((((gene_arr+i)->var1_2)/(53+120))+(((gene_arr+i)->Group_3_var)/98),0.5);
		t4 = term_4_1/term_4_2;
		term_4_3=pow((((gene_arr+i)->var1_2)/(53+120))+(((gene_arr+i)->Group_3_var)/98),2);
		term_4_4=(pow(((gene_arr+i)->var1_2)/(53+120),2)/(53+120-1))+(pow(((gene_arr+i)->Group_3_var)/98,2)/97);
		df4=term_4_3/term_4_4;
		p4_test=2*gsl_cdf_tdist_Q(fabs(t4),df4);
		if(p4_test<0.05){
			(t_test4+cnt4)->p_value=p4_test;
			strcpy((t_test4+cnt4)->passed_gene_id,(gene_arr+i)->gene_id);
			cnt4++;
		}
	}
	sort_t_test_by_p(t_test4,cnt4);
	int min4=min_get(cnt4);
	FILE *q4=fopen("Group1+2_VS_Group3.txt","w");
	for(int i=0;i<min4;i++){
		printf("GeneID:%s\tp-value:%.10lf\n",(t_test4+i)->passed_gene_id,(t_test4+i)->p_value);
		fprintf(q4,"GeneID:%s\tp-value:%.10lf\n",(t_test4+i)->passed_gene_id,(t_test4+i)->p_value);
	}
	fclose(q4);
	printf("\n\n");
	//Group2+3 VS Group1
	printf("Group2+3 VS Group1\n\n");
	struct t_test* t_test5=(struct t_test*)malloc(gene_nums*sizeof(struct t_test));
	double p5_test,term_5_1,term_5_2,t5,term_5_3,term_5_4,df5;
	int cnt5=0;
	for(int i=0;i<gene_nums;i++){
		term_5_1=(gene_arr+i)->mean2_3-(gene_arr+i)->Group_1_mean_value;
		term_5_2=pow((((gene_arr+i)->var2_3)/(120+98))+(((gene_arr+i)->Group_1_var)/53),0.5);
		t5 = term_5_1/term_5_2;
		term_5_3=pow((((gene_arr+i)->var2_3)/(120+98))+(((gene_arr+i)->Group_1_var)/53),2);
		term_5_4=(pow(((gene_arr+i)->var2_3)/(120+98),2)/(120+98-1))+(pow(((gene_arr+i)->Group_1_var)/53,2)/52);
		df5=term_5_3/term_5_4;
		p5_test=2*gsl_cdf_tdist_Q(fabs(t5),df5);
		if(p5_test<0.05){
			(t_test5+cnt5)->p_value=p5_test;
			strcpy((t_test5+cnt5)->passed_gene_id,(gene_arr+i)->gene_id);
			cnt5++;
		}
	}
	sort_t_test_by_p(t_test5,cnt5);
	int min5=min_get(cnt5);
	FILE *q5=fopen("Group2+3_VS_Group1.txt","w");
	for(int i=0;i<min5;i++){
		printf("GeneID:%s\tp-value:%.10lf\n",(t_test5+i)->passed_gene_id,(t_test5+i)->p_value);
		fprintf(q5,"GeneID:%s\tp-value:%.10lf\n",(t_test5+i)->passed_gene_id,(t_test5+i)->p_value);
	}
	fclose(q5);
	printf("\n\n");
	//Group1+3 VS Group2
	printf("Group1+3 VS Group2\n\n");
	struct t_test* t_test6=(struct t_test*)malloc(gene_nums*sizeof(struct t_test));
	double p6_test,term_6_1,term_6_2,t6,term_6_3,term_6_4,df6;
	int cnt6=0;
	for(int i=0;i<gene_nums;i++){
		term_6_1=(gene_arr+i)->mean1_3-(gene_arr+i)->Group_2_mean_value;
		term_6_2=pow((((gene_arr+i)->var1_3)/(53+98))+(((gene_arr+i)->Group_2_var)/120),0.5);
		t6 = term_6_1/term_6_2;
		term_6_3=pow((((gene_arr+i)->var1_3)/(53+98))+(((gene_arr+i)->Group_2_var)/120),2);
		term_6_4=(pow(((gene_arr+i)->var1_3)/(53+98),2)/(53+98-1))+(pow(((gene_arr+i)->Group_2_var)/120,2)/119);
		df6=term_6_3/term_6_4;
		p6_test=2*gsl_cdf_tdist_Q(fabs(t6),df6);
		if(p6_test<0.05){
			(t_test6+cnt6)->p_value=p6_test;
			strcpy((t_test6+cnt6)->passed_gene_id,(gene_arr+i)->gene_id);
			cnt6++;
		}
	}
	sort_t_test_by_p(t_test6,cnt6);
	int min6=min_get(cnt6);
	FILE *q6=fopen("Group1+3_VS_Group2.txt","w");
	for(int i=0;i<min6;i++){
		printf("GeneID:%s\tp-value:%.10lf\n",(t_test6+i)->passed_gene_id,(t_test6+i)->p_value);
		fprintf(q6,"GeneID:%s\tp-value:%.10lf\n",(t_test6+i)->passed_gene_id,(t_test6+i)->p_value);
	}
	fclose(q6);
	printf("\n\n");
	FILE *k1=fopen("Common_Gene_in_test1_and_test3.txt","w");
	printf("Common Gene in test1 and test3\n");
	fprintf(k1,"Common Gene in test1 and test3\n");
	for(int i=0;i<min1;i++){
		for(int j =0;j<min3;j++){
			if(strcmp((t_test1+i)->passed_gene_id,(t_test3+j)->passed_gene_id)==0){
				printf("%s\n",(t_test1+i)->passed_gene_id);
				fprintf(k1,"%s\n",(t_test1+i)->passed_gene_id);
			}
		}
	}
	fclose(k1);
	printf("\n\n");
	FILE *k2=fopen("Common_Gene_in_test2_and_test3.txt","w");
	printf("Common Gene in test2 and test3\n");
	fprintf(k2,"Common Gene in test2 and test3\n");
	for(int i=0;i<min2;i++){
		for(int j=0;j<min3;j++){
			if(strcmp((t_test2+i)->passed_gene_id,(t_test3+j)->passed_gene_id)==0){
				printf("%s\n",(t_test2+i)->passed_gene_id);
				fprintf(k2,"%s\n",(t_test2+i)->passed_gene_id);
			}
		}
	}
	fclose(k2);
	printf("\n\n");
	FILE *k3=fopen("Common_Gene_in_test1_and_test2.txt","w");
	printf("Common Gene in test1 and test2\n");
	fprintf(k3,"Common Gene in test1 and test2\n");
	for(int i=0;i<min1;i++){
		for(int j=0;j<min2;j++){
			if(strcmp((t_test1+i)->passed_gene_id,(t_test2+j)->passed_gene_id)==0){
				printf("%s\n",(t_test1+i)->passed_gene_id);
				fprintf(k3,"%s\n",(t_test1+i)->passed_gene_id);
			}
		}
	}
	fclose(k3);
	printf("\n\n");
	FILE *k4=fopen("Common_Gene_in_test4_and_test5.txt","w");
	printf("Common Gene in test4 and test5\n");
	fprintf(k4,"Common Gene in test4 and test5\n");
	for(int i=0;i<min4;i++){
		for(int j=0;j<min5;j++){
			if(strcmp((t_test4+i)->passed_gene_id,(t_test5+j)->passed_gene_id)==0){
				printf("%s\n",(t_test4+i)->passed_gene_id);
				fprintf(k4,"%s\n",(t_test4+i)->passed_gene_id);
			}
		}
	}
	fclose(k4);
	printf("\n\n");
	FILE *k5=fopen("Common_Gene_in_test4_and_test6.txt","w");
	printf("Common Gene in test4 and test6\n");
	fprintf(k5,"Common Gene in test4 and test6\n");
	for(int i=0;i<min4;i++){
		for(int j=0;j<min6;j++){
			if(strcmp((t_test4+i)->passed_gene_id,(t_test6+j)->passed_gene_id)==0){
				printf("%s\n",(t_test4+i)->passed_gene_id);
				fprintf(k5,"%s\n",(t_test4+i)->passed_gene_id);
			}
		}
	}
	fclose(k5);
	printf("\n\n");
	FILE *k6=fopen("Common_Gene_in_test5_and_test6.txt","w");
	printf("Common Gene in test5 and test6\n");
	fprintf(k6,"Common Gene in test5 and test6\n");
	for(int i=0;i<min5;i++){
		for(int j=0;j<min6;j++){
			if(strcmp((t_test5+i)->passed_gene_id,(t_test6+j)->passed_gene_id)==0){
				printf("%s\n",(t_test5+i)->passed_gene_id);
				fprintf(k6,"%s\n",(t_test5+i)->passed_gene_id);
			}
		}
	}
	fclose(k6);
	printf("\n\n");
	FILE *k7=fopen("Common_Gene_in_test4_and_test5_and_test6.txt","w");
	printf("Common Gene in test4 and test5 and test6\n");
	fprintf(k7,"Common Gene in test4 and test5 and test6\n");
	for(int i=0;i<min4;i++){
		for(int j=0;j<min5;j++){
			if(strcmp((t_test4+i)->passed_gene_id,(t_test5+j)->passed_gene_id)==0){
				for(int k=0;k<min6;k++){
					if(strcmp((t_test5+k)->passed_gene_id,(t_test6+k)->passed_gene_id)==0){
						printf("%s\n",(t_test4+i)->passed_gene_id);
						fprintf(k7,"%s\n",(t_test4+i)->passed_gene_id);
					}
				}
			}
			
		}
	}
	fclose(k7);
	printf("\n\n");
	FILE *k8=fopen("Common_Gene_in_test1_and_test2_and_test3.txt","w");
	printf("Common Gene in test1 and test2 and test3\n");
	fprintf(k8,"Common Gene in test1 and test2 and test3\n");
	for(int i=0;i<min1;i++){
		for(int j=0;j<min2;j++){
			if(strcmp((t_test1+i)->passed_gene_id,(t_test2+j)->passed_gene_id)==0){
				for(int k=0;k<min3;k++){
					if(strcmp((t_test2+k)->passed_gene_id,(t_test3+k)->passed_gene_id)==0){
						printf("%s\n",(t_test1+i)->passed_gene_id);
						fprintf(k7,"%s\n",(t_test1+i)->passed_gene_id);
					}
				}
			}
			
		}
	}
	fclose(k8);
    free(sample_arr);
    free(order_arr);
    free(gene_arr);
    free(t_test1);
    free(t_test2);
    free(t_test3);
    free(t_test4);
    free(t_test5);
    free(t_test6);
    return 0;
}
