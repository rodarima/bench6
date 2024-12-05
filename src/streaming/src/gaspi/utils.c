#include <mpi.h>
#include <GASPI.h>

#ifdef TAGASPI
#include <TAGASPI.h>
#endif

#include <assert.h>
#include <math.h>

#include "utils.h"
#include "common/streaming.h"

gaspi_segment_id_t RECVSEG = 0;
gaspi_segment_id_t SENDSEG = 1;

void setupGaspiInfo(const StreamingConfiguration *conf, uint64_t size, GASPInfo *info)
{
	gaspi_number_t maxNotifications;
	CHECK(gaspi_notification_num(&maxNotifications));
	if (maxNotifications < size/conf->bs) {
		fprintf(stderr, "Error: Not enough notification ids!\n");
		exit(1);
	}

	if (nranks > 1) {
		CHECK(gaspi_segment_use(RECVSEG, conf->source,
			size*sizeof(double), GASPI_GROUP_ALL,
			GASPI_BLOCK, 0));
		CHECK(gaspi_segment_use(SENDSEG, conf->array,
			size*sizeof(double), GASPI_GROUP_ALL,
			GASPI_BLOCK, 0));
	}

	CHECK(gaspi_queue_num(&info->nqueues));

#ifdef TAGASPI
	CHECK(tagaspi_queue_group_create(0, 0, info->nqueues, GASPI_QUEUE_GROUP_POLICY_CPU_RR));
#endif
}

void freeGaspiInfo(GASPInfo *info)
{
	(void) info;

#ifdef TAGASPI
	CHECK(tagaspi_queue_group_delete(0));
#endif

	if (nranks > 1) {
		CHECK(gaspi_segment_delete(RECVSEG));
		CHECK(gaspi_segment_delete(SENDSEG));
	}
}
