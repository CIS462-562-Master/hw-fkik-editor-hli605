#include "aActor.h"

#pragma warning(disable : 4018)



/****************************************************************
*
*    	    Actor functions
*
****************************************************************/

AActor::AActor()
{
	m_pInternalSkeleton = new ASkeleton();
	m_pSkeleton = m_pInternalSkeleton;

	m_BVHController = new BVHController();
	m_BVHController->setActor(this);

	m_IKController = new IKController();
	m_IKController->setActor(this);

	// code to update additional Actor data goes here
	resetGuide();

}

AActor::AActor(const AActor* actor)
{
	*this = *actor;
}

AActor& AActor::operator = (const AActor& actor)
{
	// Performs a deep copy
	if (&actor == this)
	{
		return *this;
	}
	m_pSkeleton = actor.m_pSkeleton;

	// code to update additional Actor data goes here


	return *this;
}

AActor::~AActor()
{
	delete m_IKController;
	delete m_BVHController;
	delete m_pInternalSkeleton;

}

void AActor::clear()
{
	// looks like it is clearing more times than the number of actors.  as a result, m_pSkeleton is not defined for last case.
	m_pSkeleton->clear();

	// code to update additional Actor data goes here
}

void AActor::update()
{
	if (!m_pSkeleton->getRootNode())
		return; // Nothing loaded
	else m_pSkeleton->update();

	// code to update additional Actor data goes here

}

ASkeleton* AActor::getSkeleton()
{
	return m_pSkeleton;
}

void AActor::setSkeleton(ASkeleton* pExternalSkeleton)
{
	m_pSkeleton = pExternalSkeleton;
}

void AActor::resetSkeleton()
{
	m_pSkeleton = m_pInternalSkeleton;
}

BVHController* AActor::getBVHController()
{
	return m_BVHController;
}

IKController* AActor::getIKController()
{
	return m_IKController;
}

void AActor::updateGuideJoint(vec3 guideTargetPos)
{
	if (!m_pSkeleton->getRootNode()) { return; }

	// TODO: 
	// 1.	Set the global position of the guide joint to the global position of the root joint
	// 2.	Set the y component of the guide position to 0
	// 3.	Set the global rotation of the guide joint towards the guideTarget
	AJoint* rootNode = m_pSkeleton->getRootNode();
	vec3 rootGlobalPos = m_Guide.getGlobalTranslation() + m_Guide.getGlobalRotation() * rootNode->getGlobalTranslation();
	rootGlobalPos[1] = 0;
	m_Guide.setGlobalTranslation(rootGlobalPos);

	guideTargetPos[1] = 0;
	vec3 u1 = vec3(0, 0, 1);
	vec3 u2 = (guideTargetPos - m_Guide.getGlobalTranslation()).Normalize();
	vec3 axis = u1 ^ u2;
	double angle = acos(u1 * u2);
	mat3 newRot;
	newRot.FromAxisAngle(axis, angle);
	m_Guide.setGlobalRotation(newRot);
}

void AActor::solveFootIK(float leftHeight, float rightHeight, bool rotateLeft, bool rotateRight, vec3 leftNormal, vec3 rightNormal)
{
	if (!m_pSkeleton->getRootNode()) { return; }
	AJoint* leftFoot = m_pSkeleton->getJointByID(m_IKController->mLfootID);
	AJoint* rightFoot = m_pSkeleton->getJointByID(m_IKController->mRfootID);

	// TODO: 
	// The normal and the height given are in the world space

	// 1.	Update the local translation of the root based on the left height and the right height
	vec3 rootPos = m_pSkeleton->getRootNode()->getGlobalTranslation();
	rootPos[1] += (leftHeight + rightHeight) / 2;
	m_pSkeleton->getRootNode()->setLocalTranslation(rootPos);
	
	m_pSkeleton->update();

	// 2.	Update the character with Limb-based IK 

	// Rotate Foot
	if (rotateLeft)
	{
		// Update the local orientation of the left foot based on the left normal
		ATarget m_target;
		vec3 leftFootGlobal = leftFoot->getGlobalTranslation();
		leftFootGlobal[1] = leftHeight;
		m_target.setGlobalTranslation(leftFootGlobal);

		vec3 leftNormalLocal = (leftFoot->getLocal2Global().m_rotation.Inverse() * leftNormal).Normalize();
		vec3 u = vec3(0, 1, 0);
		vec3 axis = u ^ leftNormalLocal;
		double angle = acos(u * leftNormalLocal.Normalize());
		mat3 newRot;
		newRot.FromAxisAngle(axis, angle);
		leftFoot->setLocalRotation(newRot);

		m_IKController->IKSolver_Limb(leftFoot->getID(), m_target);
	}
	if (rotateRight)
	{
		// Update the local orientation of the right foot based on the right normal
		ATarget m_target;
		vec3 rightFootGlobal = rightFoot->getGlobalTranslation();
		rightFootGlobal[1] = rightHeight;
		m_target.setGlobalTranslation(rightFootGlobal);
		

		vec3 rightNormalLocal = (rightFoot->getLocal2Global().m_rotation.Inverse() * rightNormal).Normalize();
		vec3 u = vec3(0, 1, 0);
		vec3 axis = u ^ rightNormalLocal;
		double angle = acos(u * rightNormalLocal);
		mat3 newRot;
		newRot.FromAxisAngle(axis, angle);
		rightFoot->setLocalRotation(newRot);

		m_IKController->IKSolver_Limb(rightFoot->getID(), m_target);
	}
	m_pSkeleton->update();
}
