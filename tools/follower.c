#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char** argv)
{
	if(argc != 3){
		printf("Usage: follower <file name> <timeout secs>\n");
		return -1;
	}
	struct stat s;
	size_t prev_offset = 0;
	size_t prev_ts = time(NULL);
	time_t timeout = atoi(argv[2]); // 5 minutes
	while(1){
		stat(argv[1], &s);
		size_t new_size = s.st_size - prev_offset;
		if(new_size > 0){
			FILE* ff = fopen(argv[1],"r");
			fseek(ff, prev_offset, SEEK_SET);
			char* buf = (char*)malloc(new_size);
			fread(buf, new_size, 1, ff);
			prev_offset = s.st_size;
			fwrite(buf, new_size, 1, stdout);
		}

		time_t cur_time = time(NULL);
		if(cur_time - s.st_mtime > timeout){
			break;
		}
		sleep(2);
	}
}

